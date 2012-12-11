/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
#include <mitsuba/core/sse.h>
#include <mitsuba/core/ssemath.h>
#include "../medium/materials.h"
#include "irrtree.h"
#include "bluenoise.h"

MTS_NAMESPACE_BEGIN

// Modified diffusion coefficient of Grosjean 1956 - "A high accuracy approximation for solving multiple scattering problems in infinite homogeneous media"
template <typename T> T DGrosjean(const T muA, const T muSPrime) {
    T denomterm = muA + muSPrime;
    return (2 * muA + muSPrime) / (3 * denomterm * denomterm);
}

/**
 * Computes the combined diffuse radiant exitance
 * caused by a number of dipole sources using the improved dipole model [d'Eon 2012] - "A Better Dipole".
 * This implements the R(r) given in Table 2 and weights it by the area and irradiance
 */
struct BetterDipoleQuery {
	inline BetterDipoleQuery(const Spectrum &zr,
								const Spectrum &zv,
								const Spectrum &muTr,
								const Spectrum &muA,
								const Spectrum &muSPrime,
								const Spectrum &D,
								const float &in_CE,
								const Float &in_Cphi,
								const Point &p)
		: zr(zr), zv(zv), muTr(muTr), muA(muA), muSPrime(muSPrime), m_D(D), Cphi(in_Cphi), CE(in_CE), result(0.0f), p(p) {
	}

	inline void operator()(const IrradianceSample &sample) {
		Spectrum rSqr = Spectrum((p - sample.p).lengthSquared());

		/* Distance to the real source */
		Spectrum dr = (rSqr + zr*zr).sqrt();

		/* Distance to the image point source */
		Spectrum dv = (rSqr + zv*zv).sqrt();

		Spectrum muTPrime = muA + muSPrime;

		Spectrum alpha = muSPrime / muTPrime; // single-scattering albedo

		Spectrum coeff0 = alpha * alpha; // one from Grosjean's modified diffusion, and one from the fraction of energy leaving the first scatter (at depth 1 mfp)
		Spectrum coeff1 = (CE * zr * (muTr * dr + Spectrum(1))) / (dr * dr) + Spectrum(Cphi) / m_D;
		Spectrum coeff2 = (CE * zv * (muTr * dv + Spectrum(1))) / (dv * dv) + Spectrum(Cphi) / m_D;

		/* Do not include the reduced albedo - will be canceled out later */
		Spectrum R = Spectrum(INV_FOURPI) * coeff0 *
			 (coeff1 * ((-muTr * dr).exp()) / dr
			- coeff2 * ((-muTr * dv).exp()) / dv);

		result += R * sample.E * sample.area;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	const Spectrum &zr, &zv, &muTr, &muA, &muSPrime, &m_D;
	Float Cphi;
	Float CE;
	Spectrum result;

	Point p;
};

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/*
 * Add proper exitance calculation and boundary conditions to the dipole model
 */

class BetterDipole : public Subsurface {
public:
	BetterDipole(const Properties &props)
		: Subsurface(props) {
		{
			LockGuard lock(irrOctreeMutex);
			m_octreeIndex = irrOctreeIndex++;
		}

		/* How many samples should be taken when estimating
		   the irradiance at a given point in the scene? */
		m_irrSamples = props.getInteger("irrSamples", 16);

		/* When estimating the irradiance at a given point,
		   should indirect illumination be included in the final estimate? */
		m_irrIndirect = props.getBoolean("irrIndirect", true);

		/* Multiplicative factor, which can be used to adjust the number of
		   irradiance samples */
		m_sampleMultiplier = props.getFloat("sampleMultiplier", 1.0f);

		/* Error threshold - lower means better quality */
		m_quality = props.getFloat("quality", 0.2f);

		/* Asymmetry parameter of the phase function */
		m_ready = false;
		m_octreeResID = -1;

		lookupMaterial(props, m_muS, m_muA, m_g, &m_eta);

		if (m_eta < 1)
			Log(EError, "Unsupported material configuration (intIOR/extIOR < 1)");
	}

	BetterDipole(Stream *stream, InstanceManager *manager)
	 : Subsurface(stream, manager) {
		m_muS = Spectrum(stream);
		m_muA = Spectrum(stream);
		m_g = Spectrum(stream);
		m_eta = stream->readFloat();
		m_sampleMultiplier = stream->readFloat();
		m_quality = stream->readFloat();
		m_octreeIndex = stream->readInt();
		m_irrSamples = stream->readInt();
		m_irrIndirect = stream->readBool();
		m_ready = false;
		m_octreeResID = -1;
		configure();
	}

	virtual ~BetterDipole() {
		if (m_octreeResID != -1)
			Scheduler::getInstance()->unregisterResource(m_octreeResID);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		if (m_octreeResID != -1)
			proc->bindResource(formatString("irrOctree%i", m_octreeIndex), m_octreeResID);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
		m_muS.serialize(stream);
		m_muA.serialize(stream);
		m_g.serialize(stream);
		stream->writeFloat(m_eta);
		stream->writeFloat(m_sampleMultiplier);
		stream->writeFloat(m_quality);
		stream->writeInt(m_octreeIndex);
		stream->writeInt(m_irrSamples);
		stream->writeBool(m_irrIndirect);
	}

	Spectrum LoSingle(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &dInternal) const {

		/* Sample a traveled distance within the medium */
		Float distance = -math::fastlog(1-sampler->next1D()) * m_radius;

		/* Sample a point on a light source */
		DirectSamplingRecord dRec(its.p + dInternal*distance, its.time);
		Spectrum value = scene->sampleEmitterDirect(dRec, sampler->next2D(), false);

		if (value.isZero() || !(distance > 0))
			return Spectrum(0.0f);

		/* Build a local frame based on the plane of incidence */
		Vector toLight = dRec.p - dRec.ref;
		Frame frame;
		frame.n = its.shFrame.n;
		frame.s = normalize(toLight - its.shFrame.n*dot(toLight, its.shFrame.n));
		frame.t = cross(frame.n, frame.s);

		/* Express the scattering and light source position in this frame */
		Vector V = frame.toLocal(dRec.ref - its.p);
		Vector L = frame.toLocal(dRec.p - its.p);

		if (L.z < 0)
			return Spectrum(0.0f);

		/* Find the connection path that accounts for the dielectric boundary */
		Float Vz = V.z, Lx = L.x-V.x, Lz = L.z;

		/* Starting guess proposed by Walter et al */
		Float x = (-Vz * Lx * 4) / ((m_eta * Lz - Vz) * (3 + m_eta));

		/* Perform three Newton-Raphson iterations, which is plenty */
		for (int i=0; i<3; ++i) {
			Float tmp1 = 1/(x*x + Vz*Vz), sqrtTmp1 = std::sqrt(tmp1);
			Float tmp2 = 1/((Lx-x)*(Lx-x) + Lz*Lz), sqrtTmp2 = std::sqrt(tmp2);

			Float fx  = -m_eta * x * sqrtTmp1 + (Lx - x) * sqrtTmp2;
			Float dfx = -m_eta * Vz*Vz * tmp1 * sqrtTmp1 - Lz*Lz * sqrtTmp2 * tmp2;

			x -= fx / dfx;
		}

		/* Done! */
		Vector P(V.x + x, V.y, 0);
		Point PWorld = its.p + frame.toWorld(P);

		/* Make sure that the light source is not occluded from this position */
		Vector surfaceToLight = dRec.p - PWorld;
		Float dL = surfaceToLight.length();
		surfaceToLight /= dL;

		Ray ray(PWorld, surfaceToLight, Epsilon,
			dL*(1-ShadowEpsilon), its.time);
		if (scene->rayIntersect(ray))
			return Spectrum(0.0f);

		/* Account for importance sampling wrt. transmittance */
		value *= m_radius * m_muS * (Spectrum(m_invRadius * distance) - m_muT * distance).exp();

		/* Fresnel transmittance at the new position */
		Float F = fresnelDielectricExt(surfaceToLight.z, m_eta);

		Vector interactionToSurface = P - V;
		Float dV = interactionToSurface.length();
		interactionToSurface /= dV;

		/* Evaluate the Henyey-Greenstein model */
		Float g = m_g.average(), temp = 1.0f + g*g + 2.0f * g
			* dot(interactionToSurface, -normalize(V));
		Float phase = INV_FOURPI * (1 - g*g) / (temp * std::sqrt(temp));

		Spectrum result = (1-F) * phase * value * (-m_muT * dV).exp();

		Float cosThetaL = dot(surfaceToLight, its.shFrame.n);
		Float cosThetaV = interactionToSurface.z;
		if (cosThetaL == 0 || cosThetaV == 0)
			return Spectrum(0.0f);

		/* Adjust contribution of the solution path (generalized geometric term) */
		return result * (dRec.dist*dRec.dist / (
			(dV + m_eta * dL) * (cosThetaL/cosThetaV*dV + cosThetaV/cosThetaL*m_eta*dL)));
	}

	Spectrum Lo(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		Float cosTheta = dot(its.shFrame.n, d);

		if (!m_ready || cosTheta < 0)
			return Spectrum(0.0f);

		Float F = 0, cosThetaT = 0 /* unused */;
		Vector dInternal = refract(d, its.shFrame.n, m_eta, cosThetaT, F);

#if 0
		BetterDipoleQuery query(m_zr, m_zv, m_muTr, m_muA, m_muSPrime, m_D, m_CE, m_Cphi, its.p);

		m_octree->performQuery(query);
		Spectrum result(query.getResult() * INV_PI);

		// include the normalization factor missing from Jensen et al. 2001 (see d'Eon and Irving 2011)
		if (m_eta != 1.0f)
			result *= (1.0f - F) / (1.0f - m_FdrExt);
#endif
		Spectrum result(0.0f);
		result += (1-F) * LoSingle(scene, sampler, its, dInternal);

		return result;
	}

	void configure() {
		m_muSPrime = m_muS * (Spectrum(1.0f) - m_g);
		m_muT = m_muS + m_muA;
		m_muTPrime = m_muSPrime + m_muA;
		m_invEta = 1 / m_eta;

		/* Find the smallest mean-free path over all wavelengths */
		Spectrum mfp = Spectrum(1.0f) / m_muTPrime;
		m_radius = mfp.min();
		m_invRadius = 1/m_radius;

		/* Average diffuse reflectance due to mismatched indices of refraction */
		m_FdrInt = fresnelDiffuseReflectance(1 / m_eta); // Fdr == 2 C1
		m_FdrExt = fresnelDiffuseReflectance(m_eta);

		Float _3C2 = fresnelDiffuseReflectanceSecondMoment(1 / m_eta);
		m_D = DGrosjean(m_muA, m_muSPrime);

		m_Cphi = 0.25f * (1.0f - m_FdrInt);
		m_CE  = 0.5f * (1.0f - _3C2);

		/* Dipole boundary condition distance term */
		Float A = (1 + _3C2) / (1 - m_FdrInt);

		/* Effective transport extinction coefficient */
		m_muTr = (m_muA / m_D).sqrt();

		Spectrum zb = 2.0 * A * m_D;

		/* Distance of the two dipole point sources to the surface */
		m_zr = mfp;
		m_zv = -m_zr - 2.0 * zb;
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int _samplerResID) {
		if (m_ready)
			return true;
		if (!scene->getIntegrator()->getClass()
				->derivesFrom(MTS_CLASS(SamplingIntegrator)))
			Log(EError, "The dipole subsurface scattering model requires "
				"a sampling-based surface integrator!");

		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Timer> timer = new Timer();

		AABB aabb;
		Float sa;

		ref<PositionSampleVector> points = new PositionSampleVector();
		/* It is necessary to increase the sampling resolution to
		   prevent low-frequency noise in the output */
		Float actualRadius = m_radius / std::sqrt(m_sampleMultiplier * 20);
		blueNoisePointSet(scene, m_shapes, actualRadius, points, sa, aabb, job);

		/* 2. Gather irradiance in parallel */
		const Sensor *sensor = scene->getSensor();
		ref<IrradianceSamplingProcess> proc = new IrradianceSamplingProcess(
			points, 1024, m_irrSamples, m_irrIndirect,
			sensor->getShutterOpen() + 0.5f * sensor->getShutterOpenTime(), job);

		/* Create a sampler instance for every core */
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		std::vector<SerializableObject *> samplers(sched->getCoreCount());
		for (size_t i=0; i<sched->getCoreCount(); ++i) {
			ref<Sampler> clonedSampler = sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		int samplerResID = sched->registerMultiResource(samplers);
		int integratorResID = sched->registerResource(
			const_cast<Integrator *>(scene->getIntegrator()));

		proc->bindResource("scene", sceneResID);
		proc->bindResource("integrator", integratorResID);
		proc->bindResource("sampler", samplerResID);
		scene->bindUsedResources(proc);
		m_proc = proc;
		sched->schedule(proc);
		sched->wait(proc);
		m_proc = NULL;
		for (size_t i=0; i<samplers.size(); ++i)
			samplers[i]->decRef();

		sched->unregisterResource(samplerResID);
		sched->unregisterResource(integratorResID);
		if (proc->getReturnStatus() != ParallelProcess::ESuccess)
			return false;

		Log(EDebug, "Done gathering (took %i ms), clustering ..", timer->getMilliseconds());
		timer->reset();

		std::vector<IrradianceSample> &samples = proc->getIrradianceSampleVector()->get();
		sa /= samples.size();

		for (size_t i=0; i<samples.size(); ++i)
			samples[i].area = sa;

		m_octree = new IrradianceOctree(aabb, m_quality, samples);

		Log(EDebug, "Done clustering (took %i ms).", timer->getMilliseconds());
		m_octreeResID = Scheduler::getInstance()->registerResource(m_octree);

		m_ready = true;
		return true;
	}

	void wakeup(ConfigurableObject *parent,
		std::map<std::string, SerializableObject *> &params) {
		std::string octreeName = formatString("irrOctree%i", m_octreeIndex);
		if (!m_octree.get() && params.find(octreeName) != params.end()) {
			m_octree = static_cast<IrradianceOctree *>(params[octreeName]);
			m_ready = true;
		}
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_proc);
	}

	MTS_DECLARE_CLASS()
private:
	Float m_radius, m_invRadius, m_sampleMultiplier;
	Float m_FdrInt, m_FdrExt, m_Cphi;
	Float m_CE, m_quality, m_eta, m_invEta;
	Spectrum m_muS, m_muA, m_muT;
	Spectrum m_D, m_g;
	Spectrum m_muTr, m_zr, m_zv;
	Spectrum m_muSPrime, m_muTPrime;
	ref<IrradianceOctree> m_octree;
	ref<ParallelProcess> m_proc;
	int m_octreeResID, m_octreeIndex;
	int m_irrSamples;
	bool m_irrIndirect;
	bool m_ready;
};

MTS_IMPLEMENT_CLASS_S(BetterDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(BetterDipole, "Better dipole model");
MTS_NAMESPACE_END
