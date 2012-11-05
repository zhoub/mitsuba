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
template <typename T> T DGrosjean(const T muA, const T muSPrime)
{
    T denomterm = muA + muSPrime;
    return (2 * muA + muSPrime) / (3 * denomterm * denomterm);
}

/**
 * Computes the combined diffuse radiant exitance
 * caused by a number of dipole sources
 */
struct IsotropicDipoleQuery {
	inline IsotropicDipoleQuery(const Spectrum &zr,
								const Spectrum &zv,
								const Spectrum &muTr,
								const Spectrum &muA,
								const Spectrum &muSPrime,
								const Spectrum &DCE,
								const Float &Cphi,
								const Point &p)
		: zr(zr), zv(zv), muTr(muTr), muA(muA), muSPrime(muSPrime), DCE(DCE), Cphi(Cphi), result(0.0f), p(p) {
	}

	inline void operator()(const IrradianceSample &sample) {
		Spectrum rSqr = Spectrum((p - sample.p).lengthSquared());

		/* Distance to the real source */
		Spectrum dr = (rSqr + zr*zr).sqrt();

		/* Distance to the image point source */
		Spectrum dv = (rSqr + zv*zv).sqrt();

		Spectrum muTPrime = muA + muSPrime;

		Spectrum coeff0 = 3 * muSPrime * muTPrime / (2 * muA + muSPrime);
		Spectrum coeff1 = ( DCE * zr * (muTr * dr + Spectrum(1))) / (dr * dr) + Spectrum( Cphi );
		Spectrum coeff2 = ( DCE * zv * (muTr * dv + Spectrum(1))) / (dv * dv) + Spectrum( Cphi );

		/* Do not include the reduced albedo - will be canceled out later */
		Spectrum dMo = Spectrum(INV_FOURPI) * coeff0 *
			 (coeff1 * ((-muTr * dr).exp()) / dr
			- coeff2 * ((-muTr * dv).exp()) / dv);

		result += dMo * sample.E * sample.area;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	const Spectrum &zr, &zv, &muTr, &muA, &muSPrime, &DCE;
	Float Cphi;
	Spectrum result;

	Point p;
};

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/*
 * Add proper exitance calculation and boundary conditions to the dipole model
 */

class IsotropicDipole : public Subsurface {
public:
	IsotropicDipole(const Properties &props)
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

		m_ready = false;
		m_octreeResID = -1;

		lookupMaterial(props, m_muS, m_muA, m_g, &m_eta);
	}

	IsotropicDipole(Stream *stream, InstanceManager *manager)
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

	virtual ~IsotropicDipole() {
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

	Spectrum Lo(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		if (!m_ready || dot(its.shFrame.n, d) < 0)
			return Spectrum(0.0f);
		IsotropicDipoleQuery query(m_zr, m_zv, m_muTr, m_muA, m_muSPrime, m_DCE, m_Cphi, its.p);

		m_octree->performQuery(query);
		Spectrum result(query.getResult() * INV_PI);

		// include the normalization factor missing from Jensen et al. 2001 (see d'Eon and Irving 2011)
		if (m_eta != 1.0f)
			result *= (1.0f - fresnelDielectricExt(dot(its.shFrame.n, d), m_eta)) / ( 1.0f - fresnelDiffuseReflectance( m_eta ) );

		return result;
	}

	void configure() {
		m_muSPrime = m_muS * (Spectrum(1.0f) - m_g);
		m_muTPrime = m_muSPrime + m_muA;

		/* Find the smallest mean-free path over all wavelengths */
		Spectrum mfp = Spectrum(1.0f) / m_muTPrime;
		m_radius = std::numeric_limits<Float>::max();
		for (int lambda=0; lambda<SPECTRUM_SAMPLES; lambda++)
			m_radius = std::min(m_radius, mfp[lambda]);

		/* Average diffuse reflectance due to mismatched indices of refraction */
		m_Fdr = fresnelDiffuseReflectance(1 / m_eta); // Fdr == 2 C1
		Float _3C2 = fresnelDiffuseReflectanceSecondMoment(1 / m_eta);
		Spectrum D = DGrosjean( m_muA, m_muSPrime );

		m_Cphi = 0.25f * (1.0f - m_Fdr);
		m_DCE  = D * 0.5f * (1.0f - _3C2);

		/* Dipole boundary condition distance term */
		Float A = (1 + _3C2) / (1 - m_Fdr);

		/* Effective transport extinction coefficient */
		m_muTr = (m_muA / D).sqrt();

		/* Distance of the two dipole point sources to the surface */
		m_zr = mfp / ( Spectrum(1.0f) + m_muA ); // Durian and Rudnick 1999
		m_zv = mfp * (1.0f + 4.0f/3.0f * A);
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
	Float m_radius, m_sampleMultiplier;
	Float m_Fdr, m_Cphi, m_quality, m_eta;
	Spectrum m_muS, m_muA, m_DCE, m_g;
	Spectrum m_muTr, m_zr, m_zv;
	Spectrum m_muSPrime, m_muTPrime;
	ref<IrradianceOctree> m_octree;
	ref<ParallelProcess> m_proc;
	int m_octreeResID, m_octreeIndex;
	int m_irrSamples;
	bool m_irrIndirect;
	bool m_ready;
};

MTS_IMPLEMENT_CLASS_S(IsotropicDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicDipole, "Isotropic dipole model");
MTS_NAMESPACE_END
