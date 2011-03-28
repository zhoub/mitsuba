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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Volumetric path tracer", "Average path length", EAverage);

/**
 * Volumetric path tracer, which solves the full radiative transfer
 * equation in the presence of participating media. Estimates single
 * scattering using both phase function and luminaire sampling and
 * combines the two with multiple importance sampling and the power
 * heuristic. Afterwards, the phase function sample is reused to
 * recursively estimate the multiple scattering component, which 
 * saves an intersection computation.
 * On surfaces, this integrator behaves exactly like the standard
 * MI path tracer.
 */
class VolumetricPathTracer : public MonteCarloIntegrator {
public:
	VolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) { }

	/// Unserialize from a binary data stream
	VolumetricPathTracer(Stream *stream, InstanceManager *manager)
	 : MonteCarloIntegrator(stream, manager) { }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Spectrum Li(0.0f);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;

		Spectrum pathThroughput(1.0f);

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */
			if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
				const PhaseFunction *phase = rRec.medium->getPhaseFunction();

				if (rRec.depth == m_maxDepth && m_maxDepth > 0) // No more scattering events allowed
					break;

				/* Sample the integral
				   \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/
				pathThroughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the single scattering component if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance && 
					scene->sampleAttenuatedLuminaire(mRec.p, ray.time, rRec.medium, lRec, rRec.nextSample2D())) {
					/* Evaluate the phase function */
					Spectrum phaseVal = phase->f(PhaseFunctionQueryRecord(mRec, -ray.d, -lRec.d));

					if (!phaseVal.isZero()) {
						/* Calculate prob. of having sampled that direction using 
						   phase function sampling */
						Float phasePdf = (lRec.luminaire->isIntersectable() 
								|| lRec.luminaire->isBackgroundLuminaire()) ? 
							phase->pdf(PhaseFunctionQueryRecord(mRec, -ray.d, -lRec.d)) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(lRec.pdf, phasePdf);
						Li += pathThroughput * lRec.value * phaseVal * weight;
					}
				}

				/* ==================================================================== */
				/*                         Phase function sampling                      */
				/* ==================================================================== */

				Float phasePdf;
				PhaseFunctionQueryRecord pRec(mRec, -ray.d);
				Spectrum phaseVal = phase->sample(pRec, phasePdf, rRec.sampler);
				if (phaseVal.isZero())
					break;
				phaseVal /= phasePdf;

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, pRec.wo, ray.time);
				bool hitLuminaire = false;
				if (scene->rayIntersect(ray, its)) {
					/* Intersected something - check if it was a luminaire */
					if (its.isLuminaire()) {
						lRec = LuminaireSamplingRecord(its, -ray.d);
						lRec.value = its.Le(-ray.d);
						hitLuminaire = true;
					}
				} else {
					/* No intersection found. Possibly, there is a background
					   luminaire such as an environment map? */
					if (scene->hasBackgroundLuminaire()) {
						lRec.luminaire = scene->getBackgroundLuminaire();
						lRec.d = -ray.d;
						lRec.value = lRec.luminaire->Le(ray);
						hitLuminaire = true;
					}
				}

				/* If a luminaire was hit, estimate the local illumination and
				   weight using the power heuristic */
				if (hitLuminaire && (rRec.type & RadianceQueryRecord::EDirectMediumRadiance)) {
					/* Prob. of having generated this sample using luminaire sampling */
					const Float lumPdf = scene->pdfLuminaire(mRec.p, lRec);
					Float weight = miWeight(phasePdf, lumPdf);
					Spectrum contrib = pathThroughput * lRec.value * phaseVal * weight;
					
					if (rRec.medium) {
						ray.mint = 0; ray.maxt = its.t; 
						contrib *= rRec.medium->getTransmittance(ray);
					}

					contrib += Li;
				}

				/* ==================================================================== */
				/*                         Multiple scattering                          */
				/* ==================================================================== */

				/* Stop if multiple scattering was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance)) 
					break;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;

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
					tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
					Account for this and multiply by the proper per-color-channel transmittance.
				*/
				if (rRec.medium)
					pathThroughput *= mRec.transmittance / mRec.pdfFailure;

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return 
					   attenuated radiance from a background luminaire */
					if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
						Li += pathThroughput * scene->LeBackground(ray);
					break;
				}

				/* Possibly include emitted radiance if requested */
				if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
					Li += pathThroughput * its.Le(-ray.d);

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
					Li += pathThroughput * its.LoSub(rRec.scene, -ray.d);

				if (rRec.depth == m_maxDepth && m_maxDepth > 0)
					break;

				const BSDF *bsdf = its.getBSDF(ray);
				if (!bsdf) {
					/* Pass right through the surface (there is no BSDF) */
					if (its.isMediumTransition())
						rRec.medium = its.getTargetMedium(ray.d);
					ray.setOrigin(its.p);
					scene->rayIntersect(ray, its);
					continue;
				}

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the direct illumination if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance && 
					scene->sampleAttenuatedLuminaire(its, rRec.medium, lRec, rRec.nextSample2D())) {
					/* Allocate a record for querying the BSDF */
					const BSDFQueryRecord bRec(its, its.toLocal(-lRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->fCos(bRec);

					if (!bsdfVal.isZero()) {
						/* Calculate prob. of having sampled that direction
						   using BSDF sampling */
						Float bsdfPdf = (lRec.luminaire->isIntersectable() 
								|| lRec.luminaire->isBackgroundLuminaire()) ? 
							bsdf->pdf(bRec) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(lRec.pdf, bsdfPdf);
						Li += pathThroughput * lRec.value * bsdfVal * weight;
					}
				}

				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */

				/* Sample BSDF * cos(theta) */
				BSDFQueryRecord bRec(its);
				Float bsdfPdf;
				Spectrum bsdfVal = bsdf->sampleCos(bRec, bsdfPdf, rRec.nextSample2D());
				if (bsdfVal.isZero())
					break;
				bsdfVal /= bsdfPdf;
				Intersection prevIts = its;

				/* Trace a ray in this direction */
				ray = Ray(its.p, its.toWorld(bRec.wo), ray.time);
				bool hitLuminaire = false;
				if (scene->rayIntersect(ray, its)) {
					/* Intersected something - check if it was a luminaire */
					if (its.isLuminaire()) {
						lRec = LuminaireSamplingRecord(its, -ray.d);
						lRec.value = its.Le(-ray.d);
						hitLuminaire = true;
					}
				} else {
					/* No intersection found. Possibly, there is a background
					   luminaire such as an environment map? */
					if (scene->hasBackgroundLuminaire()) {
						lRec.luminaire = scene->getBackgroundLuminaire();
						lRec.value = lRec.luminaire->Le(ray);
						lRec.d = -ray.d;
						hitLuminaire = true;
					}
				}

				/* If a luminaire was hit, estimate the local illumination and
				   weight using the power heuristic */
				if (hitLuminaire && rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) {
					/* Prob. of having generated this sample using luminaire sampling */
					const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
						scene->pdfLuminaire(prevIts.p, lRec) : 0;
					const Float weight = miWeight(bsdfPdf, lumPdf);
					Spectrum contrib = pathThroughput * lRec.value * bsdfVal * weight;
					if (rRec.medium) {
						ray.mint = 0; ray.maxt = its.t; 
						contrib *= rRec.medium->getTransmittance(ray);
					}
					Li += contrib;
				}

				/* ==================================================================== */
				/*                         Indirect illumination                        */
				/* ==================================================================== */
			
				/* Stop if indirect illumination was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) 
					break;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;

				/* Russian roulette - Possibly stop the recursion. Don't use for transmission
				   due to IOR weighting factors, which throw the heuristic off */
				if (rRec.depth >= m_rrDepth && !(bRec.sampledType & BSDF::ETransmission)) {
					/* Assuming that BSDF importance sampling is perfect,
					   'bsdfVal.max()' should equal the maximum albedo
					   over all spectral samples */
					Float approxAlbedo = std::min((Float) 0.9f, bsdfVal.max());
					if (rRec.nextSample1D() > approxAlbedo)
						break;
					else
						pathThroughput /= approxAlbedo;
				}

				pathThroughput *= bsdfVal;
				rRec.depth++;
			}
		}
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;
		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "VolumetricPathTracer[" << std::endl
			<< "  maxDepth = " << m_maxDepth << "," << std::endl
			<< "  rrDepth = " << m_rrDepth << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(VolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricPathTracer, "Volumetric path tracer");
MTS_NAMESPACE_END
