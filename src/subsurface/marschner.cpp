/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/brent.h>
#include <boost/bind.hpp>
#include "../shapes/hair.h"

MTS_NAMESPACE_BEGIN

class MarschnerModel : public Subsurface {
public:
	MarschnerModel(const Properties &props)
		: Subsurface(props) {
	}

	MarschnerModel(Stream *stream, InstanceManager *manager)
	 : Subsurface(stream, manager) {
		configure();
	}

	virtual ~MarschnerModel() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID) {
		if (!scene->getIntegrator()->getClass()->derivesFrom(MTS_CLASS(SampleIntegrator)))
			Log(EError, "Must be used with a SampleIntegrator!");
		return true;
	}

	void configure() {
		/* Precompute certain things if necessary */
	}

	/// Set the parent object
	void setParent(ConfigurableObject *parent) {
		/// Can't use derivesFrom() here for subtle linker/shared
		/// library reasons on windows
		if (parent->getClass()->getName() != "HairShape")
			Log(EError, "Can only be attached to a HairShape!");
	}

	void evalScattering(const Vector &w1, const Vector &w2, Spectrum &value) const {
		Float R = 0.6;
		value = Spectrum(R / (M_PI*M_PI));
	}

	static Float cdfErrT(Float xi, Float t) {
		// the objective function for sin(theta) sampling, described in sampleScattering below.
		return (2/M_PI * (std::acos(1-t) - std::sqrt(2*t-t*t) * (1 - t)) - xi) / std::sqrt(t);
	}

	bool sampleScattering(const Vector &w1, Vector &w2, const Point2 &sample, Spectrum &value, float &pdf) const {
		/* We want to sample the sphere according to sin(theta) where theta is the angle from the tangent.
		 * We'll do this by choosing a hemisphere and sampling the hemisphere,
		 * so the pdf to sample is (2/pi^2) sin(theta).
		 * The inverse cdf approach says
		 *    phi = 2pi xi1 and xi2 = 2pi int_0^theta (2/pi^2) sin x (sin x dx)
		 * The result is
		 *    xi2(th) = 2/pi (th - sin(th) cos(th)).
		 * This needs to be solved numerically.  A good problem transformation, arrived at empirically,
		 * is to use the independent variable t = 1 - cos(theta) and to solve
		 *    xi2(t) / sqrt(t) = xi2_given / sqrt(t)
		 * using the starting point t = xi2_given^(2/3).
		 */

		Float xi1 = sample.x, xi2 = sample.y;

		// Choose hemisphere
		bool posZ = xi2 < 0.5;
		xi2 = posZ ? 2 * xi2 : 2 - 2 * xi2;

		// initial guess
		Float t0 = pow(xi2, 2./3.);

		// solve for t
		BrentSolver brentSolver(100, 1e-6f);
		BrentSolver::Result result = brentSolver.solve(boost::bind(cdfErrT, xi2, _1), 0, 1, t0);
		SAssert(result.success);
		Float t = result.x;

		// compute theta and phi
		Float cos_theta = (posZ ? 1 : -1) * (1 - t);
		Float phi = 2*M_PI * xi1;

		// construct vector
		Float sin_theta = std::sqrt(1 - cos_theta*cos_theta);
		w2.x = cos_theta;
		w2.y = sin_theta * std::cos(phi);
		w2.z = sin_theta * std::sin(phi);

		// set value and pdf; they are constant
		evalScattering(w1, w2, value);
		pdf = evalPDF(w1, w2);

		return true;
	}


	float evalPDF(const Vector &w1, const Vector &w2) const {
		// PDF of w2 | w1 is sin(th) / pi^2.
		Float sin_theta = std::sqrt(1 - w2.x*w2.x);
		return sin_theta / (M_PI*M_PI);
	}

	// Return a point on the surface of the fiber that has the same s (along-fiber)
	// coordinate as the intersection point and is visible in the direction d.
	static Point visiblePointOnSegment(const HairShape *shape, const Intersection &its, Vector d) {

		//cerr << "geoFrame = " << its.geoFrame.toString() << endl;
		// the segment-frame offset from the segment axis to the required point
		Vector shadowPointLocal = its.geoFrame.toLocal(d);
		//cerr << "d = " << d.toString() << endl;//"; shadowPointLocal = " << shadowPointLocal.toString() << endl;
		shadowPointLocal.x = 0;
		shadowPointLocal = normalize(shadowPointLocal) * shape->getRadius();

		// the segment-coordinates position of the required point
		Point segPosn;
		int segId;
		shape->convertLocationData(its, segId, segPosn);
		segPosn.y = shadowPointLocal.y;
		segPosn.z = shadowPointLocal.z;
		//cerr << "segPosn = " << segPosn.toString() << endl;

		// the corresponding global point
		return shape->segmentToGlobal(its, segId, segPosn);
	}

	Spectrum Lo(const Scene *scene, Sampler *sampler, 
			const Intersection &its, const Vector &d, int depth) const {

		const HairShape *shape = static_cast<const HairShape *>(its.shape);

		// Result is direct lighting plus result of a recursive ray.
		Spectrum scatteredRadiance(0.0f);

		// Direct lighting: sample a luminaire point in solid angle measure, and weight by
		// the scattering function for the resulting incident direction.

		LuminaireSamplingRecord lRec;
		if (false && scene->sampleLuminaire(its.p, its.time, lRec, sampler->next2D(), false)) {

			Point shadowPoint = visiblePointOnSegment(shape, its, -lRec.d);

			Vector tempVec = shadowPoint - shape->getFirstVertex(its.color[0]);
			Vector axis = shape->getSegmentTangent(its.color[0]);
			tempVec -= axis * dot(tempVec, axis);
			//cerr << "shadowPoint " << shadowPoint.toString() << "; distance to axis = " << tempVec.toString() << endl;

			if (!scene->isOccluded(shadowPoint, lRec.sRec.p, its.time)) {
				Vector wo(its.shFrame.toLocal(-lRec.d));
				Spectrum scatFnValue;
				evalScattering(its.wi, wo, scatFnValue);
				Float cosTheta = std::sqrt(wo.y*wo.y + wo.z*wo.z);
				scatteredRadiance += cosTheta * scatFnValue * lRec.value;
//				for (int i = 0; i < depth; i++) cerr << " ";
//				cerr << "post direct: cosTheta = " << cosTheta
//						<< "; scatFnValue = " << scatFnValue.toString()
//						<< "; lRec.value = " << lRec.value.toString()
//						<< "; scatteredRadiance = " << scatteredRadiance.toString() << endl;
			}
		}

		// Recursive ray for all other lighting: generate random ray in solid angle measure,
		// weight by scattering function for that direction.
		Vector wo;
		Spectrum scatFnValue;
		Float pdf = 0;
		if (false && sampleScattering(its.wi, wo, sampler->next2D(), scatFnValue, pdf)) {

			// construct ray that can see out in the sampled direction
			Vector woWorld = its.toWorld(wo);
			Point basePoint = visiblePointOnSegment(shape, its, woWorld);

			/* Recursively gather radiance, but don't include emission */
			RadianceQueryRecord rRec(scene, sampler);
			rRec.newQuery(RadianceQueryRecord::ERadiance, NULL);
			//rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission, NULL);
			rRec.depth = depth + 1;
			Spectrum recursiveRadiance = static_cast<const SampleIntegrator *>(
					scene->getIntegrator())->Li(RayDifferential(basePoint, woWorld, its.time), rRec);

			Float cosTheta = std::sqrt(wo.y*wo.y + wo.z*wo.z);
			scatteredRadiance += cosTheta * recursiveRadiance * scatFnValue / pdf;
//			for (int i = 0; i < depth; i++) cerr << " ";
//			cerr << "post recurse: recursiveRadiance = " << recursiveRadiance.toString()
//					<< "; scatteredRadiance = " << scatteredRadiance.toString() << endl;
		}

		return scatteredRadiance;
	}

	MTS_DECLARE_CLASS()
private:
};

/* Matlab hemisphere sampling code for reference:

function [theta, phi] = sin_sample(xi1, xi2)

% Given (x1, x2) uniformly distributed in [0,1]^2, produce (theta, phi)
% distributed on the hemisphere with pdf (1/pi^2) sin(theta).

% Choose hemisphere
if xi2 < 0.5
    sign_z = 1;
    xi2 = 2*xi2;
else
    sign_z = -1;
    xi2 = 2 - 2*xi2;
end

% initial guess
t0 = xi2^(2/3);

% cdf as a function of t = (1 - cos(theta))
cdf_t = @(t) 2/pi * (acos(1-t) - sqrt(2*t-t^2) * (1 - t));

% solve for t
t = fzero(@(t) cdf_t(t) / sqrt(t) - xi2 / sqrt(t), t0);

% compute theta and phi
theta = acos(sign_z * (1 - t));
phi = 2*pi * xi1;

%% ---- testing code ----
N = 100000;
th = zeros(N,1);
ph = zeros(N,1);
for i = 1:N
  [th(i), ph(i)] = sin_sample(rand(), rand());
end

figure
hist(th, 100)
hold on
plot(linspace(0,pi), N * (pi/100) * 2 * sin(linspace(0,pi)).^2 / pi, 'r')
hold off

 */


MTS_IMPLEMENT_CLASS_S(MarschnerModel, false, Subsurface)
MTS_EXPORT_PLUGIN(MarschnerModel, "Marschner hair scattering model");
MTS_NAMESPACE_END
