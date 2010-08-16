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
		/* Number of samples to take using the luminaire sampling technique */
		m_luminaireSamples = props.getInteger("luminaireSamples", 1);
		/* Number of samples to take using the BSDF sampling technique */
		m_bsdfSamples = props.getInteger("bsdfSamples", 1);
		/* Beta factor for the power heuristic */
		m_beta = props.getFloat("beta", 2.0f);

		Assert(m_luminaireSamples + m_bsdfSamples > 0);
	}

	/// Unserialize from a binary data stream
	MIDirectIntegrator(Stream *stream, InstanceManager *manager)
	 : SampleIntegrator(stream, manager) {
		m_luminaireSamples = stream->readInt();
		m_bsdfSamples = stream->readInt();
		m_beta = stream->readFloat();
		configure();
	}

	void configure() {
		m_weightBSDF = 1 / (Float) m_bsdfSamples;
		m_weightLum = 1 / (Float) m_luminaireSamples;
		m_fracLum = m_luminaireSamples / (Float) (m_luminaireSamples + m_bsdfSamples);
		m_fracBSDF = m_bsdfSamples / (Float) (m_luminaireSamples + m_bsdfSamples);
	}

	void configureSampler(Sampler *sampler) {
		sampler->request2DArray(m_luminaireSamples);
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
		if (rRec.rayIntersect(ray)) {
			ray.mint = 0; ray.maxt = rRec.its.t;
		}

		if (!its.isValid()) {
			/* If no intersection could be found, possibly return 
			   attenuated radiance from a background luminaire */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				return scene->LeBackground(ray);
			else
				return Spectrum(0.0f);
		}

		const BSDF *bsdf = its.getBSDF(ray);

		/* Possibly include emitted radiance if requested */
		if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface integrator if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(scene, -ray.d);

		/* Leave here if direct illumination was not requested */
		if (!(rRec.type & RadianceQueryRecord::EDirectRadiance))
			return Li * rRec.attenuation;

		/* ==================================================================== */
		/*                          Luminaire sampling                          */
		/* ==================================================================== */
		Point2 *sampleArray;
		int numLuminaireSamples = m_luminaireSamples,
			numBSDFSamples = m_bsdfSamples;
		Float fracLum = m_fracLum, fracBSDF = m_fracBSDF,
		      weightLum = m_weightLum, weightBSDF = m_weightBSDF;

		if (rRec.depth > 1) {
			/* This integrator is used by another integrator.
			   Be less accurate as this sample will not directly be observed. */
			numBSDFSamples = numLuminaireSamples = 1;
			fracLum = fracBSDF = .5f;
			weightLum = weightBSDF = 1.0f;
		}

		if (numLuminaireSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numLuminaireSamples);
		} else {
			sample = rRec.nextSample2D();
			sampleArray = &sample;
		}

		for (int i=0; i<numLuminaireSamples; ++i) {
			/* Estimate the direct illumination if this is requested */
			if (scene->sampleLuminaire(its, lRec, sampleArray[i])) {
				/* Allocate a record for querying the BSDF */
				const BSDFQueryRecord bRec(rRec, its, its.toLocal(-lRec.d));

				/* Evaluate BSDF * cos(theta) */
				const Spectrum bsdfVal = bsdf->fCos(bRec);

				if (!bsdfVal.isBlack()) {
					/* Calculate prob. of having sampled that direction
						using BSDF sampling */
					Float bsdfPdf = (lRec.luminaire->isIntersectable() 
							|| lRec.luminaire->isBackgroundLuminaire()) ? 
						bsdf->pdf(bRec) : 0;

					/* Weight using the power heuristic */
					const Float weight = miWeight(lRec.pdf * fracLum, bsdfPdf * fracBSDF) * weightLum;
					Li += lRec.Le * bsdfVal * weight;
				}
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		} else {
			sample = rRec.nextSample2D();
			sampleArray = &sample;
			Assert(sampleArray[0] == sample);
		}

		for (int i=0; i<numBSDFSamples; ++i) {
			/* Sample BSDF * cos(theta) */
			BSDFQueryRecord bRec(rRec, its, sampleArray[i]);
			Float bsdfPdf;
			Spectrum bsdfVal = bsdf->sampleCos(bRec, bsdfPdf);
			if (bsdfVal.isBlack())
				break;
			bsdfVal /= bsdfPdf;

			/* Trace a ray in this direction */
			Ray bsdfRay(its.p, its.toWorld(bRec.wo));
			bool hitLuminaire = false;
			if (scene->rayIntersect(bsdfRay, bsdfIts)) {
				bsdfRay.mint = 0; bsdfRay.maxt = bsdfIts.t; 
				/* Intersected something - check if it was a luminaire */
				if (bsdfIts.isLuminaire()) {
					lRec = LuminaireSamplingRecord(bsdfIts, -bsdfRay.d);
					hitLuminaire = true;
				}
			} else {
				/* No intersection found. Possibly, there is a background
					luminaire such as an environment map? */
				if (scene->hasBackgroundLuminaire()) {
					lRec.luminaire = scene->getBackgroundLuminaire();
					lRec.d = -bsdfRay.d;
					hitLuminaire = true;
				}
			}

			/* If a luminaire was hit, estimate the local illumination and
				sample weight using the power heuristic */
			if (hitLuminaire && (rRec.type & RadianceQueryRecord::EDirectRadiance)) {
				lRec.Le = lRec.luminaire->Le(lRec);
				Float lumPdf = scene->pdfLuminaire(its, lRec);
				if (bRec.sampledType & BSDF::EDelta)
					lumPdf = 0;
				const Float weight = miWeight(bsdfPdf * fracBSDF, lumPdf * fracLum) * weightBSDF;
				Li += lRec.Le * bsdfVal * weight;
			}
		}

		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA = std::pow(pdfA, m_beta);
		pdfB = std::pow(pdfB, m_beta);
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SampleIntegrator::serialize(stream, manager);
		stream->writeInt(m_luminaireSamples);
		stream->writeInt(m_bsdfSamples);
		stream->writeFloat(m_beta);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIDirectIntegrator[" << std::endl
			<< "  luminaireSamples = " << m_luminaireSamples << "," << std::endl
			<< "  bsdfSamples = " << m_bsdfSamples << "," << std::endl
			<< "  beta = " << m_beta << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_luminaireSamples, m_bsdfSamples;
	Float m_fracBSDF, m_fracLum;
	Float m_weightBSDF, m_weightLum;
	Float m_beta;
};

MTS_IMPLEMENT_CLASS_S(MIDirectIntegrator, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(MIDirectIntegrator, "Direct-only integrator");
MTS_NAMESPACE_END
