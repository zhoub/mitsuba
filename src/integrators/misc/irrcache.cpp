/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/plugin.h>
#include "irrcache_proc.h"

MTS_NAMESPACE_BEGIN

/**
 * Irradiance caching integrator - forwards all radiance computations 
 * to an arbitrary sampling-based sub-integrator - with one exception:
 * whenever a Lambertian surface is intersected, an internal irradiance
 * cache is queried for the indirect illumination at the surface position in 
 * question. If this query is successful, the sub-integrator is only 
 * used to compute the remaining types of radiance (direct, in-scatter, 
 * emission) and their sum is returned afterwards.
 * When a query is unsuccessful, a new data point is generated by a final
 * gathering step.
 *
 * The generality of this implementation allows it to be used in conjunction
 * with photon mapping (the most likely application) as well as all other 
 * sampling-based integrators in Mitsuba. Several optimizations are used to
 * improve the achieved interpolation quality, namely irradiance gradients 
 * [Ward et al.], neighbor clamping [Krivanek et al.], a screen-space 
 * clamping metric and an improved error function [Tabellion et al.].
 * By default, this integrator also performs a distributed overture pass before
 * rendering, which is recommended to avoid artifacts resulting from the
 * addition of samples as rendering proceeds.
 */
class IrradianceCacheIntegrator : public SampleIntegrator {
	friend class OvertureThread;
public:
	IrradianceCacheIntegrator(const Properties &props) : SampleIntegrator(props) {
		/* Elevational resolution of the stratified final gather hemisphere.
		   The azimuthal resolution is three times this value. Default: 
		   14x(3*14)=588 samples */
		m_resolution = props.getInteger("resolution", 14);
		/* If set to true, the irradiance cache will be filled by a
		   parallel overture pass before the main rendering process starts.
		   This is strongly recommended. */
		m_overture = props.getBoolean("overture", true);
		/* Quality setting (\kappa in the [Tabellion et al.] paper).
		   A value of 1 should be adequate in most cases. */
		m_quality = props.getFloat("quality", 1.0f);
		/* Multiplicative factor for the quality parameter following an
		   overture pass. This can be used to interpolate amongst more
		   samples, creating a visually smoother result. Must be
		   1 or less.  */
		m_qualityAdjustment = props.getFloat("qualityAdjustment", .5f);
		/* If set to true, sample locations will be visually highlighted */
		m_debug = props.getBoolean("debug", false);
		/* Should irradiance gradients be used? Generally, this will
		   significantly improve the interpolation quality.*/
		m_gradients = props.getBoolean("gradients", true);
		/* Should neighbor clamping [Krivanek et al.] be used? This 
		   propagates geometry information amongst close-by samples 
		   and generally leads to better sample placement. */
		m_clampNeighbor = props.getBoolean("clampNeighbor", true);
		/* If set to true, the influence region of samples will be clamped
		   using the screen-space metric by [Tabellion et al.]? 
		   Turning this off may lead to excessive sample placement. */
		m_clampScreen = props.getBoolean("clampScreen", true);
		/* Minimum influence region of an irradiance sample (relative to scene size, in [0,1]) */
		m_influenceMin = props.getFloat("influenceMin", 0.005f);
		/* Maximum influence region of an irradiance sample (default=64*min) */
		m_influenceMax = props.getFloat("influenceMax", 64*m_influenceMin);
		/* If set to false, direct illumination will be suppressed - 
		   useful for checking the interpolation quality */
		m_direct = props.getBoolean("direct", true);
			
		Assert(m_influenceMax > m_influenceMin);
		Assert(m_influenceMax > 0 && m_influenceMax < 1);
		Assert(m_influenceMin > 0 && m_influenceMin < 1);
		Assert(m_qualityAdjustment > 0 && m_qualityAdjustment <= 1);
	}

	IrradianceCacheIntegrator(Stream *stream, InstanceManager *manager) 
	 : SampleIntegrator(stream, manager) {
		m_irrCache = static_cast<IrradianceCache *>(manager->getInstance(stream));
		m_subIntegrator = static_cast<SampleIntegrator *>(manager->getInstance(stream));
		m_resolution = stream->readInt();
		m_influenceMin = stream->readFloat();
		m_influenceMax = stream->readFloat();
		m_quality = stream->readFloat();
		m_qualityAdjustment = stream->readFloat();
		m_clampScreen = stream->readBool();
		m_clampNeighbor = stream->readBool();
		m_overture = stream->readBool();
		m_gradients = stream->readBool();
		m_debug = stream->readBool();
		m_direct = stream->readBool();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SampleIntegrator::serialize(stream, manager);
		manager->serialize(stream, m_irrCache.get());
		manager->serialize(stream, m_subIntegrator.get());
		stream->writeInt(m_resolution);
		stream->writeFloat(m_influenceMin);
		stream->writeFloat(m_influenceMax);
		stream->writeFloat(m_quality);
		stream->writeFloat(m_qualityAdjustment);
		stream->writeBool(m_clampScreen);
		stream->writeBool(m_clampNeighbor);
		stream->writeBool(m_overture);
		stream->writeBool(m_gradients);
		stream->writeBool(m_debug);
		stream->writeBool(m_direct);
	}

	void configureSampler(Sampler *sampler) {
		m_subIntegrator->configureSampler(sampler);
	}
	
	void bindUsedResources(ParallelProcess *proc) const {
		m_subIntegrator->bindUsedResources(proc);
	}

	void wakeup(std::map<std::string, SerializableObject *> &params) {
		m_subIntegrator->wakeup(params);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
			if (!cClass->derivesFrom(MTS_CLASS(SampleIntegrator)))
				Log(EError, "The sub-integrator must be derived from the class SampleIntegrator");
			m_subIntegrator = static_cast<SampleIntegrator *>(child);
		} else {
			Integrator::addChild(name, child);
		}
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID) {
		if (!SampleIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID))
			return false;

		if (m_subIntegrator == NULL)
			Log(EError, "No sub-integrator was specified!");

		if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID))
			return false;

		ref<Scheduler> sched = Scheduler::getInstance();
		m_irrCache = new IrradianceCache(scene->getAABB());
		m_irrCache->clampNeighbor(m_clampNeighbor);
		m_irrCache->clampScreen(m_clampScreen);
		m_irrCache->clampInfluence(m_influenceMin, m_influenceMax);
		m_irrCache->useGradients(m_gradients);
		m_irrCache->setQuality(m_quality);

		std::string irrCacheStatus;
		if (m_overture)
			irrCacheStatus += "overture, ";
		if (m_debug)
			irrCacheStatus += "debug, ";
		if (m_gradients)
			irrCacheStatus += "gradients, ";
		if (m_clampNeighbor)
			irrCacheStatus += "clampNeighbor, ";
		if (m_clampScreen)
			irrCacheStatus += "clampScreen, ";
		irrCacheStatus += formatString("clampWorld(%.3f,%.3f)", m_influenceMin, m_influenceMax);

		Log(EDebug, "Irradiance cache status : %s", irrCacheStatus.c_str());
		Log(EDebug, "  - Gather resolution   : %ix%i = %i samples", m_resolution, 3*m_resolution, 3*m_resolution*m_resolution);
		Log(EDebug, "  - Quality setting     : %.2f (adjustment: %.2f)", m_quality, m_qualityAdjustment);

		if (m_overture) {
			int subIntegratorResID = sched->registerResource(m_subIntegrator);
			ref<OvertureProcess> proc = new OvertureProcess(job, m_resolution, m_gradients, 
				m_clampNeighbor, m_clampScreen, m_influenceMin, m_influenceMax, m_quality);
			m_proc = proc;
			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("subIntegrator", subIntegratorResID);
			bindUsedResources(proc);
			sched->schedule(proc);
			sched->unregisterResource(subIntegratorResID);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
				Log(EWarn, "The overture pass did not complete sucessfully!");
				return false;
			}

			ref<const IrradianceRecordVector> vec = proc->getSamples();
			Log(EDebug, "Overture pass generated %i irradiance samples", vec->size());
			for (size_t i=0; i<vec->size(); ++i)
				m_irrCache->insert(new IrradianceCache::Record((*vec)[i]));

			m_irrCache->setQuality(m_quality * m_qualityAdjustment);
		}
		return true;
	}

	void cancel() {
		if (m_proc) {
			Scheduler::getInstance()->cancel(m_proc);
		} else {
			SampleIntegrator::cancel();
			m_subIntegrator->cancel();
		}
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		Intersection &its = rRec.its;
		if (!m_direct && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance))
			rRec.type ^= RadianceQueryRecord::EDirectSurfaceRadiance;

		if (rRec.rayIntersect(ray)) {
			const BSDF *bsdf = its.getBSDF(ray);

			if (bsdf->getType() == BSDF::EDiffuseReflection && 
					(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
				Spectrum E;
				if (!m_irrCache->get(its, E)) {
					handleMiss(ray, rRec, E);

					if (m_debug)
						E.fromLinearRGB(1e3, 0, 0);
				}

				rRec.type ^= RadianceQueryRecord::EIndirectSurfaceRadiance;

				return E * bsdf->getDiffuseReflectance(its) * INV_PI + 
					m_subIntegrator->Li(ray, rRec);
			}
		}
		return m_subIntegrator->Li(ray, rRec);
	}

	void handleMiss(const RayDifferential &ray, const RadianceQueryRecord &rRec,
			Spectrum &E) const {
		/* Handle an irradiance cache miss */
		HemisphereSampler *hs = m_hemisphereSampler.get();
		Sampler *sampler = m_sampleGenerator.get();
		RadianceQueryRecord rRec2;
		if (hs == NULL) {
			Properties props("independent");
			props.setInteger("sampleCount", m_resolution * 3 * m_resolution);
			sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), props));
			hs = new HemisphereSampler(m_resolution, 3 * m_resolution);
			m_hemisphereSampler.set(hs);
			m_sampleGenerator.set(sampler);
		}

		/* Generate stratified cosine-weighted samples and compute
		   rotational + translational gradients */
		hs->generateDirections(rRec.its, sampler);
		sampler->generate();

		for (unsigned int j=0; j<hs->getM(); j++) {
			for (unsigned int k=0; k<hs->getN(); k++) {
				HemisphereSampler::SampleEntry &entry = (*hs)(j, k);
					entry.dist = std::numeric_limits<Float>::infinity();
				rRec2.recursiveQuery(rRec, 
					RadianceQueryRecord::ERadianceNoEmission | RadianceQueryRecord::EDistance);
				rRec2.extra = 1;
				rRec2.sampler = sampler;
				entry.L = m_subIntegrator->Li(RayDifferential(rRec.its.p, entry.d, ray.time), rRec2);
				entry.dist = rRec2.dist;
				sampler->advance();
			}
		}

		hs->process(rRec.its);
		m_irrCache->put(ray, rRec.its, *hs);
		E = hs->getIrradiance();
	}

	Spectrum E(const Scene *scene, const Point &p, const Normal &n, Float time,
			const Medium *medium, Sampler *sampler, int nSamples, bool handleIndirect) const {
		Spectrum EDir(0.0f), EIndir(0.0f);
		RadianceQueryRecord rRec(scene, sampler);
		LuminaireSamplingRecord lRec;

		/* Direct illumination */
		for (int i=0; i<nSamples; i++) {
			if (scene->sampleAttenuatedLuminaire(p, time, medium, lRec, rRec.nextSample2D())) {
				Float dp = dot(lRec.d, n);
				if (dp < 0) 
					EDir -= lRec.value * dp;
			}
		}

		rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission, medium);
		rRec.its.p = p;
		rRec.its.geoFrame = rRec.its.shFrame = Frame(n);

		if (handleIndirect) {
			if (!m_irrCache->get(rRec.its, EIndir)) 
				handleMiss(RayDifferential(), rRec, EIndir);
		}

		return (EDir / (Float) nSamples) + EIndir;
	}

	const Integrator *getSubIntegrator() const {
		return m_subIntegrator.get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "IrradianceCacheIntegrator[" << std::endl
			<< "  subIntegrator = " << indent(m_subIntegrator->toString()) << "," << std::endl
			<< "  resolution = " << m_resolution << "," << std::endl
			<< "  irrCache = " << indent(m_irrCache->toString()) << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	mutable ThreadLocal<HemisphereSampler> m_hemisphereSampler;
	mutable ThreadLocal<Sampler> m_sampleGenerator;
	mutable ref<IrradianceCache> m_irrCache;
	ref<SampleIntegrator> m_subIntegrator;
	ref<ParallelProcess> m_proc;
	int m_resolution;
	Float m_influenceMin, m_influenceMax;
	Float m_quality, m_qualityAdjustment;
	bool m_clampScreen, m_clampNeighbor;
	bool m_overture, m_gradients, m_debug, m_direct;
};

MTS_IMPLEMENT_CLASS_S(IrradianceCacheIntegrator, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(IrradianceCacheIntegrator, "Irradiance cache");
MTS_NAMESPACE_END
