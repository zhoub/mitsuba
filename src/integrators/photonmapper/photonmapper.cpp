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

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/gatherproc.h>
#include "bre.h"

MTS_NAMESPACE_BEGIN

/**
* Parallel photon gathering, SSE-accelerated lookups, RGBE encoding
* Numbers for caustic/volume are ignored when no specular objects or participating
* media are present.
*/
class PhotonMapIntegrator : public SampleIntegrator {
public:
	PhotonMapIntegrator(const Properties &props) : SampleIntegrator(props) {
		/* Number of luminaire samples for direct illumination */
		m_directSamples = props.getInteger("directSamples", 16);
		/* Number of BSDF samples when intersecting a glossy material */
		m_glossySamples = props.getInteger("glossySamples", 32);
		/* Depth to start using russian roulette when tracing photons */
		m_rrDepth = props.getInteger("rrDepth", 10);
		/* Depth cutoff when tracing photons */
		m_maxDepth = props.getInteger("maxDepth", 40);
		/* Depth cutoff when recursively tracing specular materials */
		m_maxSpecularDepth = props.getInteger("maxSpecularDepth", 6);
		/* Granularity of photon tracing work units (in shot particles, 0 => decide automatically) */
		m_granularity = props.getInteger("granularity", 0);
		/* Number of photons to collect for the global photon map */
		m_globalPhotons = props.getSize("globalPhotons", 200000);
		/* Number of photons to collect for the caustic photon map */
		m_causticPhotons = props.getSize("causticPhotons", 200000);
		/* Number of photons to collect for the volumetric photon map */
		m_volumePhotons = props.getSize("volumePhotons", 200000);
		/* Radius of lookups in the global photon map (relative to the scene size) */
		m_globalLookupRadiusRel = props.getFloat("globalLookupRadius", 0.05f);
		/* Radius of lookups in the caustic photon map (relative to the scene size) */
		m_causticLookupRadiusRel = props.getFloat("causticLookupRadius", 0.0125f);
		/* Minimum amount of photons to consider a volumetric photon map lookup valid */
		m_globalLookupSize = props.getInteger("globalLookupSize", 120);
		/* Maximum number of results for caustic photon map lookups */
		m_causticLookupSize = props.getInteger("causticLookupSize", 120);
		/* Approximate number of volume photons to be used in a lookup */
		m_volumeLookupSize = props.getInteger("volumeLookupSize", 120);
		/* Should photon gathering steps exclusively run on the local machine? */
		m_gatherLocally = props.getBoolean("gatherLocally", true);
	}

	/// Unserialize from a binary data stream
	PhotonMapIntegrator(Stream *stream, InstanceManager *manager)
	 : SampleIntegrator(stream, manager) {
		m_directSamples = stream->readInt();
		m_glossySamples = stream->readInt();
		m_maxSpecularDepth = stream->readInt();
		m_globalPhotons = stream->readSize();
		m_causticPhotons = stream->readSize();
		m_volumePhotons = stream->readSize();
		m_globalLookupRadius = stream->readFloat();
		m_causticLookupRadius = stream->readFloat();
		m_globalLookupSize = stream->readInt();
		m_causticLookupSize = stream->readInt();
		m_volumeLookupSize = stream->readInt();
		m_gatherLocally = stream->readBool();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SampleIntegrator::serialize(stream, manager);
		stream->writeInt(m_directSamples);
		stream->writeInt(m_glossySamples);
		stream->writeInt(m_maxSpecularDepth);
		stream->writeSize(m_globalPhotons);
		stream->writeSize(m_causticPhotons);
		stream->writeSize(m_volumePhotons);
		stream->writeFloat(m_globalLookupRadius);
		stream->writeFloat(m_causticLookupRadius);
		stream->writeInt(m_globalLookupSize);
		stream->writeInt(m_causticLookupSize);
		stream->writeInt(m_volumeLookupSize);
		stream->writeBool(m_gatherLocally);
	}

	/// Configure the sampler for a specified amount of direct illumination samples
	void configureSampler(Sampler *sampler) {
		if (m_directSamples > 1)
			sampler->request2DArray(m_directSamples);
		sampler->request2DArray(m_glossySamples);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, 
			int sceneResID, int cameraResID, int samplerResID) {
		SampleIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);
		/* Create a deterministic sampler for the photon gathering step */
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("halton")));
		int qmcSamplerID = sched->registerResource(sampler);

		/* Don't create a caustic photon map if the scene does not contain specular materials */
		const std::vector<Shape *> &shapes = scene->getShapes();
		bool foundSpecular = false;
		for (size_t i=0; i<shapes.size(); ++i) {
			const BSDF *bsdf = shapes[i]->getBSDF();
			if (bsdf && bsdf->getType() & BSDF::EDelta) {
				foundSpecular = true;
				break;
			}
		}
		if (!foundSpecular)
			m_causticPhotons = 0;

		/* Don't create a volumetric photon map if there are no participating media */
		const std::set<Medium *> &media = scene->getMedia();
		if (media.size() == 0)
			m_volumePhotons = 0;

		for (std::set<Medium *>::const_iterator it = media.begin(); it != media.end(); ++it) {
			if (!(*it)->isHomogeneous())
				Log(EError, "Inhomogeneous media are currently not supported by the photon mapper!");
		}


		if (m_globalPhotonMap.get() == NULL && m_globalPhotons > 0) {
			/* Adapt to scene extents */
			m_globalLookupRadius = m_globalLookupRadiusRel * scene->getBSphere().radius;

			/* Generate the global photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::ESurfacePhotons, m_globalPhotons,
				m_granularity, m_maxDepth, m_rrDepth, m_gatherLocally, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

			m_globalPhotonMap = proc->getPhotonMap();
			m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
			m_globalPhotonMap->balance();
			m_globalPhotonMapID = sched->registerResource(m_globalPhotonMap);
		}

		if (m_causticPhotonMap.get() == NULL && m_causticPhotons > 0) {
			/* Adapt to scene extents */
			m_causticLookupRadius = m_causticLookupRadiusRel * scene->getBSphere().radius;

			/* Generate the caustic photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::ECausticPhotons, m_causticPhotons,
				m_granularity, 2, m_rrDepth, m_gatherLocally, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;
	
			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;
			
			Log(EDebug, "Caustic photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

			m_causticPhotonMap = proc->getPhotonMap();
			m_causticPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
			m_causticPhotonMap->balance();
			m_causticPhotonMapID = sched->registerResource(m_causticPhotonMap);
		}

		if (m_volumePhotonMap.get() == NULL && m_volumePhotons > 0) {
			/* Generate the volume photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::EVolumePhotons, m_volumePhotons,
				m_granularity, m_maxDepth, m_rrDepth, m_gatherLocally, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			Log(EDebug, "Volume photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

			ref<PhotonMap> volumePhotonMap = proc->getPhotonMap();
			volumePhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
			volumePhotonMap->balance();
	
			m_bre = new BeamRadianceEstimator(volumePhotonMap, m_volumeLookupSize);
			m_breID = sched->registerResource(m_bre);
		}

		sched->unregisterResource(qmcSamplerID);
		m_parentIntegrator = static_cast<SampleIntegrator *>(getParent());
		return true;
	}

	/// Specify globally shared resources
	void bindUsedResources(ParallelProcess *proc) const {
		if (m_globalPhotonMap.get())
			proc->bindResource("globalPhotonMap", m_globalPhotonMapID);
		if (m_causticPhotonMap.get())
			proc->bindResource("causticPhotonMap", m_causticPhotonMapID);
		if (m_bre.get())
			proc->bindResource("bre", m_breID);
	}

	/// Connect to globally shared resources
	void wakeup(std::map<std::string, SerializableObject *> &params) {
		if (!m_globalPhotonMap.get() && params.find("globalPhotonMap") != params.end())
			m_globalPhotonMap = static_cast<PhotonMap *>(params["globalPhotonMap"]);
		if (!m_causticPhotonMap.get() && params.find("causticPhotonMap") != params.end())
			m_causticPhotonMap = static_cast<PhotonMap *>(params["causticPhotonMap"]);
		if (!m_bre.get() && params.find("bre") != params.end())
			m_bre = static_cast<BeamRadianceEstimator *>(params["bre"]);

		if (getParent() != NULL && getParent()->getClass()->derivesFrom(MTS_CLASS(SampleIntegrator)))
			m_parentIntegrator = static_cast<SampleIntegrator *>(getParent());
		else
			m_parentIntegrator = this;

	}

	void cancel() {
		SampleIntegrator::cancel();
		if (m_proc)
			Scheduler::getInstance()->cancel(m_proc);
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		Spectrum LiSurf(0.0f), LiMedium(0.0f), transmittance(1.0f);
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
	
		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);

		if (rRec.medium) {
			Ray mediumRaySegment(ray, 0, its.t);
			transmittance = rRec.medium->getTransmittance(mediumRaySegment);
			mediumRaySegment.mint = ray.mint;
			if (rRec.type & RadianceQueryRecord::EVolumeRadiance)
				LiMedium = m_bre->query(mediumRaySegment, rRec.medium);
		}

		if (!its.isValid()) {
			/* If no intersection could be found, possibly return 
			   attenuated radiance from a background luminaire */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				LiSurf += rRec.scene->LeBackground(ray);
			return LiSurf * transmittance + LiMedium;
		}

		/* Possibly include emitted radiance if requested */
		if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)) 
			LiSurf += its.Le(-ray.d);

		/* Include radiance from a subsurface integrator if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			LiSurf += its.LoSub(rRec.scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);
				
		if (bsdf == NULL) {
			if (rRec.depth+1 < m_maxSpecularDepth) {
				RadianceQueryRecord rRec2;
				rRec2.recursiveQuery(rRec);
	
				if (its.isMediumTransition())
					rRec2.medium = its.getTargetMedium(ray.d);

				LiSurf += m_parentIntegrator->Li(RayDifferential(its.p, ray.d, ray.time), rRec2);
			}
			return LiSurf * transmittance + LiMedium;
		}

		int bsdfType = bsdf->getType();

		Point2 *sampleArray, sample;
		int numDirectSamples = (rRec.depth == 1) ? m_directSamples : 1;
		/* When this integrator is used recursively by another integrator,
			Be less accurate as this sample will not be directly observed. */
		if (numDirectSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numDirectSamples);
		} else {
			sample = rRec.nextSample2D();
			sampleArray = &sample;
		}

		/* Estimate the direct illumination if this is requested */
		if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) {
			Float weight = 1 / (Float) numDirectSamples;

			for (int i=0; i<numDirectSamples; ++i) {
				if (rRec.scene->sampleAttenuatedLuminaire(its, rRec.medium, lRec, sampleArray[i])) {
					/* Allocate a record for querying the BSDF */
					const BSDFQueryRecord bRec(its, its.toLocal(-lRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);

					LiSurf += lRec.value * bsdfVal * weight;
				}
			}
		}

		if (bsdfType == BSDF::EDiffuseReflection) {
			/* Hit a diffuse material - do a direct photon map visualization. */
			if (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)
				LiSurf += m_globalPhotonMap->estimateIrradianceFiltered(its.p, 
					its.shFrame.n, m_globalLookupRadius, m_globalLookupSize)
						* bsdf->getDiffuseReflectance(its) * INV_PI;
			if (rRec.type & RadianceQueryRecord::ECausticRadiance && m_causticPhotonMap.get())
				LiSurf += m_causticPhotonMap->estimateIrradianceFiltered(its.p,
					its.shFrame.n, m_causticLookupRadius, m_causticLookupSize)
							* bsdf->getDiffuseReflectance(its) * INV_PI;
		} else if ((bsdfType & BSDF::EDelta) != 0
				&& (bsdfType & ~BSDF::EDelta) == 0 && !rRec.extra) {
			RadianceQueryRecord rRec2;
			RayDifferential recursiveRay;
			/* Ideal specular material -> recursive ray tracing.
			   Deliberately risk exponential ray growth by spawning
			   several child rays. The impact on the final image is huge
			   and well worth the extra computation. */
			if (rRec.depth+1 < m_maxSpecularDepth) {
				int compCount = bsdf->getComponentCount();
				for (int i=0; i<compCount; i++) {
					/* Sample the BSDF and recurse */
					BSDFQueryRecord bRec(its);
					bRec.component = i;
					Spectrum bsdfVal = bsdf->sample(bRec, Point2(0.0f));
					if (bsdfVal.isZero())
						continue;

					rRec2.recursiveQuery(rRec, RadianceQueryRecord::ERadiance);
					recursiveRay = Ray(its.p, its.toWorld(bRec.wo), ray.time);
					LiSurf += m_parentIntegrator->Li(recursiveRay, rRec2) * bsdfVal;
				}
			}
		} else if (rRec.depth == 1 && (bsdf->getType() & BSDF::EGlossy)) {
			/* Hit a glossy material - MC integration over the hemisphere
			   using BSDF importance sampling */
			sampleArray = rRec.sampler->next2DArray(m_glossySamples);
			RadianceQueryRecord rRec2;
			RayDifferential recursiveRay;
			Float weight = 1 / (Float) m_glossySamples;

			for (int i=0; i<m_glossySamples; ++i) {
				BSDFQueryRecord bRec(its);
				bRec.sampler = rRec.sampler;
				Spectrum bsdfVal = bsdf->sample(bRec, sampleArray[i]);

				rRec2.recursiveQuery(rRec, RadianceQueryRecord::ERadianceNoEmission);
				recursiveRay = Ray(its.p, its.toWorld(bRec.wo), ray.time);
				LiSurf += m_parentIntegrator->Li(recursiveRay, rRec2) * bsdfVal * weight;
			}
		} else {
			LiSurf += m_globalPhotonMap->estimateRadianceFiltered(its,
				m_globalLookupRadius, m_globalLookupSize);
		}

		return LiSurf * transmittance + LiMedium;
	}
	
	std::string toString() const {
		std::ostringstream oss;
		oss << "PhotonMapIntegrator[" << std::endl
			<< "  directSamples = " << m_directSamples << "," << std::endl
			<< "  glossySamples = " << m_glossySamples << "," << std::endl
			<< "  globalPhotons = " << m_globalPhotons << "," << std::endl
			<< "  causticPhotons = " << m_causticPhotons << "," << std::endl
			<< "  volumePhotons = " << m_volumePhotons << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<PhotonMap> m_globalPhotonMap;
	ref<PhotonMap> m_causticPhotonMap;
	ref<PhotonMap> m_volumePhotonMap;
	ref<ParallelProcess> m_proc;
	ref<BeamRadianceEstimator> m_bre;
	SampleIntegrator *m_parentIntegrator;
	int m_globalPhotonMapID, m_causticPhotonMapID, m_breID;
	size_t m_globalPhotons, m_causticPhotons, m_volumePhotons;
	int m_globalLookupSize, m_causticLookupSize, m_volumeLookupSize;
	Float m_globalLookupRadiusRel, m_globalLookupRadius;
	Float m_causticLookupRadiusRel, m_causticLookupRadius;
	int m_granularity;
	int m_directSamples, m_glossySamples;
	int m_rrDepth;
	int m_maxDepth, m_maxSpecularDepth;
	bool m_gatherLocally;
};

MTS_IMPLEMENT_CLASS_S(PhotonMapIntegrator, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(PhotonMapIntegrator, "Photon map integrator");
MTS_NAMESPACE_END
