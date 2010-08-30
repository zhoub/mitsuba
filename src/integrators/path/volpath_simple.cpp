#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Volumetric path tracer, which solves the full radiative transfer
 * equation in the presence of participating media. Simplified version
 * without multiple importance sampling - this version can be much
 * faster than the multiple importance sampling version when when 
 * rendering heterogeneous participating media using the 
 * [Coleman et al.] sampling technique, since fewer attenuation
 * evaluations will be required.
 */
class SimpleVolumetricPathTracer : public MonteCarloIntegrator {
public:
	SimpleVolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) {
	}

	/// Unserialize from a binary data stream
	SimpleVolumetricPathTracer(Stream *stream, InstanceManager *manager)
	 : MonteCarloIntegrator(stream, manager) {
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Intersection prevIts;
		Spectrum Li(0.0f);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		Spectrum pathThroughput(1.0f);
		bool computeIntersection = false;

		while (rRec.depth < m_maxDepth || m_maxDepth < 0) {
			if (computeIntersection)
				scene->rayIntersect(ray, its);

			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */
			if (scene->sampleDistance(ray, its.t, mRec, rRec.sampler)) {
				const PhaseFunction *phase = mRec.medium->getPhaseFunction();
				Vector wo, wi = -ray.d;

				/* Sample the integral
				   \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/
				pathThroughput *= mRec.sigmaS * mRec.attenuation * mRec.miWeight / mRec.pdf;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the single scattering component if this is requested */
				if (rRec.type & RadianceQueryRecord::EInscatteredDirectRadiance && 
					scene->sampleLuminaireAttenuated(mRec.p, lRec, rRec.nextSample2D())) {
					Li += pathThroughput * lRec.Le * phase->f(mRec, -ray.d, -lRec.d);
				}

				/* ==================================================================== */
				/*                         Phase function sampling                      */
				/* ==================================================================== */

				PhaseFunction::ESampledType sampledType;
				Spectrum phaseVal = phase->sample(mRec, wi, wo, sampledType, rRec.nextSample2D());
				if (phaseVal.max() == 0)
					break;
				prevIts = its;

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, wo);
				computeIntersection = true;

				/* ==================================================================== */
				/*                         Multiple scattering                          */
				/* ==================================================================== */

				/* Set the recursive query type */
				if (!(rRec.type & RadianceQueryRecord::EInscatteredIndirectRadiance)) {
					/* Stop if multiple scattering was not requested (except: sampled a delta phase function
					   - look for emitted radiance only) */
					if (sampledType == PhaseFunction::EDelta) 
						rRec.type = RadianceQueryRecord::EEmittedRadiance;
					else
						break;
				} else {
					if (sampledType != PhaseFunction::EDelta || !(rRec.type & RadianceQueryRecord::EInscatteredDirectRadiance)) {
						/* Emitted radiance is only included in the recursive query if:
						  - the sampled phase function component had a Dirac delta distribution AND
						  - the current query asks for single scattering
						*/
						rRec.type = RadianceQueryRecord::ERadianceNoEmission;
					} else {
						rRec.type = RadianceQueryRecord::ERadiance;
					}
				}

				/* Russian roulette - Possibly stop the recursion */
				if (rRec.depth >= m_rrDepth) {
					if (rRec.nextSample1D() > mRec.albedo)
						break;
					else
						pathThroughput /= mRec.albedo;
				}

				pathThroughput *= phaseVal;
				rRec.depth++;
			} else {
				/* Sample 
					tau(x, y) (Surface integral). This happens with probability mRec.pdf
					Divide this out and multiply with the proper per-color-channel attenuation.
				*/
				pathThroughput *= mRec.attenuation * mRec.miWeight / mRec.pdf;

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return 
					   attenuated radiance from a background luminaire */
					if (rRec.type & RadianceQueryRecord::EEmittedRadiance) 
						Li += pathThroughput * scene->LeBackground(ray);
					break;
				}

				const BSDF *bsdf = its.getBSDF(ray);

				/* Possibly include emitted radiance if requested */
				if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
					Li += pathThroughput * its.Le(-ray.d);

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
					Li += pathThroughput * its.LoSub(rRec.scene, -ray.d);

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the direct illumination if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectRadiance && 
					scene->sampleLuminaireAttenuated(its, lRec, rRec.nextSample2D())) {
					/* Allocate a record for querying the BSDF */
					const BSDFQueryRecord bRec(rRec, its, its.toLocal(-lRec.d));

					Li += pathThroughput * lRec.Le * bsdf->fCos(bRec);
				}

				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */

				/* Sample BSDF * cos(theta) */
				BSDFQueryRecord bRec(rRec, its, rRec.nextSample2D());
				Spectrum bsdfVal = bsdf->sampleCos(bRec);
				if (bsdfVal.isBlack())
					break;
				prevIts = its;

				/* Trace a ray in this direction */
				ray = Ray(its.p, its.toWorld(bRec.wo));
				computeIntersection = true;

				/* ==================================================================== */
				/*                         Indirect illumination                        */
				/* ==================================================================== */
			
				/* Stop if indirect illumination was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectRadiance)) {
					/* Stop if indirect illumination was not requested (except: sampled a delta BSDF 
					   - look for emitted radiance only) */
					if (bRec.sampledType & BSDF::EDelta) 
						rRec.type = RadianceQueryRecord::EEmittedRadiance;
					else
						break;
				} else {
					if (!(bRec.sampledType & BSDF::EDelta) || !(rRec.type & RadianceQueryRecord::EDirectRadiance)) {
						/* Emitted radiance is only included in the recursive query if:
						  - the sampled BSDF component had a Dirac delta distribution AND
						  - the current query asks for direct illumination
						*/
						rRec.type = RadianceQueryRecord::ERadianceNoEmission;
					} else {
						rRec.type = RadianceQueryRecord::ERadiance;
					}
				}

				/* Russian roulette - Possibly stop the recursion */
				if (rRec.depth >= m_rrDepth) {
					/* Assuming that BSDF importance sampling is perfect,
					   the following should equal the maximum albedo
					   over all spectral samples */
					Float approxAlbedo = std::min((Float) 0.9, bsdfVal.max());
					if (rRec.nextSample1D() > approxAlbedo)
						break;
					else
						pathThroughput /= approxAlbedo;
				}

				pathThroughput *= bsdfVal;
				rRec.depth++;
			}
		}
		return Li;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SimpleVolumetricPathTracer[" << std::endl
			<< "  maxDepth = " << m_maxDepth << "," << std::endl
			<< "  rrDepth = " << m_rrDepth << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_beta;
};

MTS_IMPLEMENT_CLASS_S(SimpleVolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SimpleVolumetricPathTracer, "Simple volumetric path tracer");
MTS_NAMESPACE_END
