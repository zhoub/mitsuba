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

MTS_NAMESPACE_BEGIN

/**
 * Direct-only integrator using multiple importance sampling and
 * the power heuristic. Takes a user-specifiable amount of luminaire
 * and BSDF samples By setting one of the strategies to zero, this 
 * class can effectively be turned into a luminaire sampling or
 * BSDF sampling-based integrator. Ignores participating media.
 */
class MIDirectIntegrator : public SampleIntegrator {
public:
	MIDirectIntegrator(const Properties &props) : SampleIntegrator(props) {
		/* Number of shading samples -- this parameter is a shorthand notation
		   to set both 'luminaireSamples' and 'bsdfSamples' at the same time*/
		size_t shadingSamples = props.getSize("shadingSamples", 1);

		/* Number of samples to take using the luminaire sampling technique */
		m_luminaireSamples = props.getSize("luminaireSamples", shadingSamples);
		/* Number of samples to take using the BSDF sampling technique */
		m_bsdfSamples = props.getSize("bsdfSamples", shadingSamples);

		Assert(m_luminaireSamples >= 0 && m_bsdfSamples >= 0 &&
			m_luminaireSamples + m_bsdfSamples > 0);
	}

	/// Unserialize from a binary data stream
	MIDirectIntegrator(Stream *stream, InstanceManager *manager)
	 : SampleIntegrator(stream, manager) {
		m_luminaireSamples = stream->readSize();
		m_bsdfSamples = stream->readSize();
		configure();
	}

	void configure() {
		size_t sum = m_luminaireSamples + m_bsdfSamples;
		m_weightBSDF = 1 / (Float) m_bsdfSamples;
		m_weightLum = 1 / (Float) m_luminaireSamples;
		m_fracBSDF = m_bsdfSamples / (Float) sum;
		m_fracLum = m_luminaireSamples / (Float) sum;
	}

	void configureSampler(Sampler *sampler) {
		if (m_luminaireSamples > 1)
			sampler->request2DArray(m_luminaireSamples);
		if (m_bsdfSamples > 1)
			sampler->request2DArray(m_bsdfSamples);
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its, bsdfIts;
		RayDifferential ray(r);
		LuminaireSamplingRecord lRec;
		Spectrum Li(0.0f);
		Point2 sample;

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		if (!rRec.rayIntersect(ray)) {
			/* If no intersection could be found, possibly return 
			   radiance from a background luminaire */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				return scene->LeBackground(ray);
			else
				return Spectrum(0.0f);
		}

		/* Possibly include emitted radiance if requested */
		if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface integrator if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);

		if (EXPECT_NOT_TAKEN(!bsdf)) {
			/* The direct illumination integrator doesn't support
			   surfaces without a BSDF (e.g. medium transitions)
			   -- give up. */
			return Li;
		}

		/* Leave here if direct illumination was not requested */
		if (!(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance))
			return Li;

		/* ==================================================================== */
		/*                          Luminaire sampling                          */
		/* ==================================================================== */
		Point2 *sampleArray;
		size_t numLuminaireSamples = m_luminaireSamples,
			   numBSDFSamples = m_bsdfSamples;
		Float fracLum = m_fracLum, fracBSDF = m_fracBSDF,
		      weightLum = m_weightLum, weightBSDF = m_weightBSDF;

		if (rRec.depth > 1) {
			/* This integrator is used recursively by another integrator.
			   Be less accurate as this sample will not directly be observed. */
			numBSDFSamples = numLuminaireSamples = 1;
			fracLum = fracBSDF = .5f;
			weightLum = weightBSDF = 1.0f;
		}

		if (numLuminaireSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numLuminaireSamples);
		} else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		for (size_t i=0; i<numLuminaireSamples; ++i) {
			/* Estimate the direct illumination if this is requested */
			if (scene->sampleLuminaire(its.p, ray.time, lRec, sampleArray[i])) {
				/* Allocate a record for querying the BSDF */
				BSDFQueryRecord bRec(its, its.toLocal(-lRec.d));

				/* Evaluate BSDF * cos(theta) */
				const Spectrum bsdfVal = bsdf->eval(bRec);

				if (!bsdfVal.isZero()) {
					/* Calculate prob. of having sampled that direction
						using BSDF sampling */
					Float bsdfPdf = (lRec.luminaire->isIntersectable() 
							|| lRec.luminaire->isBackgroundLuminaire()) ? 
						bsdf->pdf(bRec) : 0;

					/* Weight using the power heuristic */
					const Float weight = miWeight(lRec.pdf * fracLum, 
							bsdfPdf * fracBSDF) * weightLum;
					Li += lRec.value * bsdfVal * weight;
				}
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		} else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		for (size_t i=0; i<numBSDFSamples; ++i) {
			/* Sample BSDF * cos(theta) */
			BSDFQueryRecord bRec(its, rRec.sampler, ERadiance);
			Float bsdfPdf;
			Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
			if (bsdfVal.isZero())
				continue;

			/* Trace a ray in this direction */
			Ray bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);
			if (scene->rayIntersect(bsdfRay, bsdfIts)) {
				/* Intersected something - check if it was a luminaire */
				if (bsdfIts.isLuminaire()) {
					lRec = LuminaireSamplingRecord(bsdfIts, -bsdfRay.d);
					lRec.value = bsdfIts.Le(-bsdfRay.d);
				} else {
					continue;
				}
			} else {
				/* No intersection found. Possibly, there is a background
				   luminaire such as an environment map? */
				if (scene->hasBackgroundLuminaire()) {
					lRec.luminaire = scene->getBackgroundLuminaire();
					lRec.d = -bsdfRay.d;
					lRec.value = lRec.luminaire->Le(bsdfRay);
				} else {
					continue;
				}
			}

			const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ? 
				scene->pdfLuminaire(its.p, lRec) : 0;
	
			const Float weight = miWeight(bsdfPdf * fracBSDF, 
				lumPdf * fracLum) * weightBSDF;
			Li += lRec.value * bsdfVal * weight;
		}

		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA; pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SampleIntegrator::serialize(stream, manager);
		stream->writeSize(m_luminaireSamples);
		stream->writeSize(m_bsdfSamples);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIDirectIntegrator[" << endl
			<< "  luminaireSamples = " << m_luminaireSamples << "," << endl
			<< "  bsdfSamples = " << m_bsdfSamples << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	size_t m_luminaireSamples;
	size_t m_bsdfSamples;
	Float m_fracBSDF, m_fracLum;
	Float m_weightBSDF, m_weightLum;
};

MTS_IMPLEMENT_CLASS_S(MIDirectIntegrator, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(MIDirectIntegrator, "Direct illumination integrator");
MTS_NAMESPACE_END
