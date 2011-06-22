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
		value = Spectrum(1/(2*M_PI));
	}

	bool sampleScattering(const Vector &w1, Vector &w2, Spectrum &value, float &pdf) const {
		return false;
	}

	float evalPDF(const Vector &w1, const Vector &w2) const {
		return 1/(2*M_PI);
	}

	// Return a point on the surface of the fiber that has the same s (along-fiber)
	// coordinate as the intersection point and is visible in the direction d.
	static Point visiblePointOnSegment(const HairShape *shape, const Intersection &its, Vector d) {

		//cerr << "geoFrame = " << its.geoFrame.toString() << endl;
		// the segment-frame offset from the segment axis to the required point
		Vector shadowPointLocal = its.geoFrame.toLocal(-d);
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
		if (scene->sampleLuminaire(its.p, its.time, lRec, sampler->next2D(), false)) {

			Point shadowPoint = visiblePointOnSegment(shape, its, lRec.d);

			Vector tempVec = shadowPoint - shape->getFirstVertex(its.color[0]);
			Vector axis = shape->getSegmentTangent(its.color[0]);
			tempVec -= axis * dot(tempVec, axis);
			//cerr << "shadowPoint " << shadowPoint.toString() << "; distance to axis = " << tempVec.toString() << endl;

			if (!scene->isOccluded(shadowPoint, lRec.sRec.p, its.time)) {
				Vector wo(its.shFrame.toLocal(lRec.d));
				Spectrum scatFnValue;
				evalScattering(its.wi, wo, scatFnValue);
				Float cosTheta = std::sqrt(wo.y*wo.y + wo.z*wo.z);
				scatteredRadiance += cosTheta * scatFnValue * lRec.value / lRec.pdf;
			}
		}

		// Recursive ray for all other lighting: generate random ray in solid angle measure,
		// weight by scattering function for that direction.
		Vector wo;
		Spectrum scatFnValue;
		Float pdf = 0;
		if (sampleScattering(its.wi, wo, scatFnValue, pdf)) {

			/* Recursively gather radiance, but don't include emission */
			RadianceQueryRecord rRec(scene, sampler);
			rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission, NULL);
			rRec.depth = depth + 1;

			if (sampleScattering(its.wi, wo, scatFnValue, pdf)) {

				Point basePoint = visiblePointOnSegment(shape, its, wo);

				Spectrum recursiveRadiance = static_cast<const SampleIntegrator *>(
						scene->getIntegrator())->Li(RayDifferential(basePoint, wo, its.time), rRec);

				Float cosTheta = std::sqrt(wo.y*wo.y + wo.z*wo.z);
				scatteredRadiance += cosTheta * recursiveRadiance * scatFnValue / pdf;
			}
		}

		return scatteredRadiance;
	}

	MTS_DECLARE_CLASS()
private:
};

MTS_IMPLEMENT_CLASS_S(MarschnerModel, false, Subsurface)
MTS_EXPORT_PLUGIN(MarschnerModel, "Marschner hair scattering model");
MTS_NAMESPACE_END
