#include <mitsuba/render/vpl.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

size_t generateVPLs(const Scene *scene, size_t offset, size_t count, int maxDepth, 
		std::deque<VPL> &vpls) {
	ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
		createObject(Sampler::m_theClass, Properties("halton")));
	EmissionRecord eRec;
	eRec.type = EmissionRecord::EPreview;
	Ray ray;
	Intersection its;
	Spectrum weight, bsdfVal;
	int depth;

	if (maxDepth <= 1)
		return 0;

	const Frame stdFrame(Vector(1,0,0), Vector(0,1,0), Vector(0,0,1));

	while (vpls.size() < count) {
		sampler->setSampleIndex(++offset);
		Point2 areaSample = sampler->next2D(),
		       dirSample  = sampler->next2D();

		/* Sample an emitted particle */
		scene->sampleEmissionArea(eRec, areaSample);
		weight = eRec.P / eRec.pdfArea;
		VPL lumVPL(ELuminaireVPL, weight);
		lumVPL.its.p = eRec.sRec.p;
		lumVPL.its.shFrame = eRec.luminaire->getType() == Luminaire::EDeltaPosition 
			? stdFrame : Frame(eRec.sRec.n);
		lumVPL.luminaire = eRec.luminaire;
		vpls.push_back(lumVPL);

		weight *= scene->sampleEmissionDirection(eRec, dirSample);
		Float cosTheta = eRec.sRec.n.isZero() ? (Float) 1 : absDot(eRec.sRec.n, eRec.d);
		weight *= cosTheta / eRec.pdfDir;
		ray = Ray(eRec.sRec.p, eRec.d);

		depth = 2;
		while (!weight.isBlack() && depth < maxDepth) {
			if (!scene->rayIntersect(ray, its))
				break;

			const BSDF *bsdf = its.shape->getBSDF();
			BSDFQueryRecord bRec(its, sampler->next2D());
			bRec.quantity = EImportance;
			bsdfVal = bsdf->sampleCos(bRec);
			if (bsdfVal.isBlack())
				break;

			/* Assuming that BSDF importance sampling is perfect,
				the following should equal the maximum albedo
				over all spectral samples */
			Float approxAlbedo = std::min((Float) 1, bsdfVal.max());
			if (sampler->next1D() > approxAlbedo)
				break;
			else
				weight /= approxAlbedo;

			VPL vpl(ESurfaceVPL, weight);
			vpl.its = its;
			vpls.push_back(vpl);
	
			weight *= bsdfVal;
		
			Vector wi = -ray.d, wo = its.toWorld(bRec.wo);
			ray = Ray(its.p, wo);

			/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
			Float wiDotGeoN = dot(its.geoFrame.n, wi),
				  woDotGeoN = dot(its.geoFrame.n, wo);
			if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 || 
				woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;

			/* Adjoint BSDF for shading normals -- [Veach, p. 155] */
			weight *= std::abs(
				(Frame::cosTheta(bRec.wi) * woDotGeoN)/
				(Frame::cosTheta(bRec.wo) * wiDotGeoN));

			++depth;
		}
	}
	return offset;
}

const char *toString(EVPLType type) {
	switch (type) {
		case ELuminaireVPL: return "luminaireVPL";
		case ESurfaceVPL: return "surfaceVPL";
		default:
			SLog(EError, "Unknown VPL type!");
			return NULL;
	}
}

std::string VPL::toString() const {
	std::ostringstream oss;
	oss << "VPL[" << endl
		<< "  type = " << mitsuba::toString(type) << "," << endl
		<< "  P = " << P.toString() << "," << endl;
	if (type == ELuminaireVPL) {
		oss << "  p = " << its.p.toString() << "," << endl;
		oss << "  luminaire = " << indent(luminaire->toString()) << endl;
	} else {
		oss << "  its = " << indent(its.toString()) << endl;
	}
	oss << "]" << endl;
	return oss.str();
}

MTS_NAMESPACE_END
