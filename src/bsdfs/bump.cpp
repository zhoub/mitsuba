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
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{bump}{Bump map}
 *
 * \parameters{
 *     \parameter{\Unnamed}{\Texture}{
 *       The luminance of this texture specifies the amount of 
 *       displacement. The implementation ignores any constant
 *       offset---only changes in the luminance matter.
 *     }
 *     \parameter{\Unnamed}{\BSDF}{A BSDF model that should
 *     be affected by the bump map}
 * }
 * \renderings{
 *     \rendering{Bump map based on tileable diagonal lines}{bsdf_bump_1}
 *     \rendering{An irregular bump map}{bsdf_bump_2}
 * }
 *
 * Bump mapping \cite{Blinn1978Simulation} is a simple technique for cheaply
 * adding surface detail to a rendering. This is done by perturbing the
 * shading coordinate frame based on a displacement height field provided
 * as a texture. This method can lend objects a highly realistic and detailed 
 * appearance (e.g. wrinkled or covered by scratches and other imperfections)
 * without requiring any changes to the input geometry.
 *
 * The implementation in Mitsuba uses the common approach of ignoring
 * the usually negligible texture-space derivative of the base mesh 
 * surface normal. As side effect of this decision, it is invariant
 * to constant offsets in the height field texture---only variations in
 * its luminance cause changes to the shading frame.
 *
 * Note that the magnitude of the height field variations influences 
 * the strength of the displacement. If desired, the \pluginref{scale}
 * texture plugin can be used to magnify or reduce the effect of a 
 * bump map texture.
 * \begin{xml}[caption=A rough metal model with a scaled image-based bump map]
 * <bsdf type="bump">
 *     <!-- The bump map is applied to a rough metal BRDF -->
 *     <bsdf type="roughconductor"/>
 *
 *     <texture type="scale">
 *         <!-- The scale of the displacement gets multiplied by 10x -->
 *         <float name="scale" value="10"/>
 *
 *         <texture type="bitmap">
 *             <string name="filename" value="bumpmap.png"/>
 *         </texture>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
class BumpMap : public BSDF {
public:
	BumpMap(const Properties &props) : BSDF(props) { }

	BumpMap(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_displacement = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");
		if (!m_displacement)
			Log(EError, "A displacement texture must be specified");

		m_components.clear();
		for (int i=0; i<m_nested->getComponentCount(); ++i) 
			m_components.push_back(m_nested->getType(i) | ESpatiallyVarying | EAnisotropic);

		m_usesRayDifferentials = true;

		BSDF::configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_nested.get());
		manager->serialize(stream, m_displacement.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (m_nested != NULL)
				Log(EError, "Only a single nested BSDF can be added!");
			m_nested = static_cast<BSDF *>(child);
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (m_displacement != NULL)
				Log(EError, "Only a single displacement texture can be specified!");
			m_displacement = static_cast<Texture *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void perturbIntersection(const Intersection &its, Intersection &target) const {
		/* Determine the step size for the finite difference computation */
		Float du = 0.5f * std::abs(its.dudx) + std::abs(its.dudy),
			  dv = 0.5f * std::abs(its.dvdx) + std::abs(its.dvdy);
		if (du == 0) du = Epsilon;
		if (dv == 0) dv = Epsilon;

		/* Compute the U and V displacementment derivatives */
		Float displacement, displacementU, displacementV;
		target = its;
		displacement = m_displacement->getValue(target).getLuminance();
		target.p = its.p + its.dpdu * du;
		target.uv = its.uv + Point2(du, 0);
		displacementU = m_displacement->getValue(target).getLuminance();
		target.p = its.p + its.dpdv * dv;
		target.uv = its.uv + Point2(0, dv);
		displacementV = m_displacement->getValue(target).getLuminance();
		target.p = its.p;
		target.uv = its.uv;
		Float dDisplaceDu = (displacementU - displacement) / du;
		Float dDisplaceDv = (displacementV - displacement) / dv;

		/* Build a perturbed frame -- ignores the usually 
		   negligible normal derivative term */
		Vector dpdu = its.dpdu + its.shFrame.n * (
				dDisplaceDu - dot(its.shFrame.n, its.dpdu));
		Vector dpdv = its.dpdv + its.shFrame.n * (
				dDisplaceDv - dot(its.shFrame.n, its.dpdv));

		dpdu = normalize(dpdu);
		dpdv = normalize(dpdv - dpdu * dot(dpdv, dpdu));

		target.shFrame.s = dpdu;
		target.shFrame.t = dpdv;
		target.shFrame = Frame(Normal(cross(dpdv, dpdu)));

		if (dot(target.shFrame.n, target.geoFrame.n) < 0)
			target.shFrame.n *= -1;
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		const Intersection& its = bRec.its;
		Intersection perturbed;
		perturbIntersection(its, perturbed);

		BSDFQueryRecord perturbedQuery(perturbed,
			perturbed.toLocal(its.toWorld(bRec.wi)),
			perturbed.toLocal(its.toWorld(bRec.wo)), bRec.quantity);
		if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
			return Spectrum(0.0f);
		perturbedQuery.sampler = bRec.sampler;
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		return m_nested->eval(perturbedQuery, measure);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		const Intersection& its = bRec.its;
		Intersection perturbed;
		perturbIntersection(its, perturbed);
		
		BSDFQueryRecord perturbedQuery(perturbed,
			perturbed.toLocal(its.toWorld(bRec.wi)),
			perturbed.toLocal(its.toWorld(bRec.wo)), bRec.quantity);
		if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
			return 0;
		perturbedQuery.quantity = bRec.quantity;
		perturbedQuery.sampler = bRec.sampler;
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		return m_nested->pdf(perturbedQuery, measure);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		const Intersection& its = bRec.its;
		Intersection perturbed;
		perturbIntersection(its, perturbed);

		BSDFQueryRecord perturbedQuery(perturbed, bRec.sampler, bRec.quantity);
		perturbedQuery.wi = perturbed.toLocal(its.toWorld(bRec.wi));
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		Spectrum result = m_nested->sample(perturbedQuery, pdf, sample);

		if (!result.isZero()) {
			bRec.sampledComponent = perturbedQuery.sampledComponent;
			bRec.sampledType = perturbedQuery.sampledType;
			bRec.wo = its.toLocal(perturbed.toWorld(perturbedQuery.wo));
			if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
				return Spectrum(0.0f);
		}

		return result;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		const Intersection& its = bRec.its;
		Intersection perturbed;
		perturbIntersection(its, perturbed);

		BSDFQueryRecord perturbedQuery(perturbed, bRec.sampler, bRec.quantity);
		perturbedQuery.wi = perturbed.toLocal(its.toWorld(bRec.wi));
		perturbedQuery.sampler = bRec.sampler;
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		Spectrum result = m_nested->sample(perturbedQuery, sample);
		if (!result.isZero()) {
			bRec.sampledComponent = perturbedQuery.sampledComponent;
			bRec.sampledType = perturbedQuery.sampledType;
			bRec.wo = its.toLocal(perturbed.toWorld(perturbedQuery.wo));
			if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
				return Spectrum(0.0f);
		}
		return result;
	}

	Shader *createShader(Renderer *renderer) const;

	std::string toString() const {
		std::ostringstream oss;
		oss << "BumpMap[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  displacement = " << indent(m_displacement->toString()) << endl
			<< "  nested = " << indent(m_nested->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	ref<Texture> m_displacement;
	ref<BSDF> m_nested;
};

// ================ Hardware shader implementation ================ 

/**
 * This is a quite approximate version of the bump map model -- it likely
 * won't match the reference exactly, but it should be good enough for
 * preview purposes
 */
class BumpMapShader : public Shader {
public:
	BumpMapShader(Renderer *renderer, const BSDF *nested, const Texture *displacement) 
		: Shader(renderer, EBSDFShader), m_nested(nested), m_displacement(displacement) {
		m_nestedShader = renderer->registerShaderForResource(m_nested.get());
		m_displacementShader = renderer->registerShaderForResource(m_displacement.get());
	}

	bool isComplete() const {
		return m_nestedShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nested.get());
		renderer->unregisterShaderForResource(m_displacement.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedShader.get());
		deps.push_back(m_displacementShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    float du = abs(dFdx(uv.x)), dv = abs(dFdx(uv.y));" << endl
			<< "    if (du == 0.0) du = 0.001;" << endl
			<< "    if (dv == 0.0) dv = 0.001;" << endl
			<< "    float displacement = " << depNames[1] << "(uv)[0];" << endl
			<< "    float displacementU = " << depNames[1] << "(uv + vec2(du, 0.0))[0];" << endl
			<< "    float displacementV = " << depNames[1] << "(uv + vec2(0.0, dv))[0];" << endl
			<< "    float dfdu = (displacementU - displacement)/du;" << endl
			<< "    float dfdv = (displacementV - displacement)/dv;" << endl
			<< "    vec3 dpdu = normalize(vec3(1.0, 0.0, dfdu));" << endl
			<< "    vec3 dpdv = vec3(0.0, 1.0, dfdv);" << endl
			<< "    dpdv = normalize(dpdv - dot(dpdu, dpdv)*dpdu);" << endl
			<< "    vec3 n = cross(dpdu, dpdv);" << endl
			<< "    wi = vec3(dot(wi, dpdu), dot(wi, dpdv), dot(wi, n));" << endl
			<< "    wo = vec3(dot(wo, dpdu), dot(wo, dpdv), dot(wo, n));" << endl
			<< "    return " << depNames[0] << "(uv, wi, wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const BSDF> m_nested;
	ref<const Texture> m_displacement;
	ref<Shader> m_nestedShader;
	ref<Shader> m_displacementShader;
};

Shader *BumpMap::createShader(Renderer *renderer) const { 
	return new BumpMapShader(renderer, m_nested.get(), m_displacement.get());
}

MTS_IMPLEMENT_CLASS(BumpMapShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(BumpMap, false, BSDF)
MTS_EXPORT_PLUGIN(BumpMap, "Smooth dielectric coating");
MTS_NAMESPACE_END
