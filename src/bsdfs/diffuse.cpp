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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{diffuse}{Smooth diffuse material}
 *
 * \parameters{
 *     \lastparameter{reflectance}{\Spectrum\Or\Texture}{
 *       Specifies the diffuse reflectance / albedo of the material \linebreak(Default: 0.5)
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Homogeneous reflectance, see \lstref{diffuse-uniform}}{bsdf_diffuse_plain}
 *     \rendering{Textured reflectance, see \lstref{diffuse-textured}}{bsdf_diffuse_textured}
 * }
 *
 * The smooth diffuse material (often referred to as ``Lambertian'' materials)
 * represents an ideally diffuse material with a user-specified amount of 
 * reflectance. Any received illumination is scattered so that the surface 
 * looks the same independently of the direction of observation.
 *
 * Apart from a  homogeneous reflectance value, the plugin can also accept 
 * a nested or referenced texture map to be used as the source of reflectance 
 * information, which is then mapped onto the shape based on its UV
 * parameterization. When no parameters are specified, the model uses the default 
 * of 50% reflectance.
 *
 * Note that this material is one-sided---that is, observed from the 
 * back side, it will be completely black. If this is undesirable, 
 * consider using the \pluginref{twosided} BRDF adapter plugin.
 *
 * \begin{xml}[caption=Reflectance specified as an sRGB color, label=lst:diffuse-uniform]
 * <bsdf type="diffuse">
 *     <srgb name="reflectance" value="#6d7185"/>
 * </bsdf>
 * \end{xml}
 *
 * \begin{xml}[caption=Diffuse material with a texture map, label=lst:diffuse-textured]
 * <bsdf type="diffuse">
 *     <texture type="bitmap" name="reflectance">
 *         <string name="filename" value="wood.jpg"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
class SmoothDiffuse : public BSDF {
public:
	SmoothDiffuse(const Properties &props) 
		: BSDF(props) {
		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance" 
				: "diffuseReflectance", Spectrum(.5f)));
		m_components.push_back(EDiffuseReflection | EFrontSide);
		m_usesRayDifferentials = false;
	}

	SmoothDiffuse(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_components.push_back(EDiffuseReflection | EFrontSide);
		m_usesRayDifferentials = m_reflectance->usesRayDifferentials();
	}

	virtual ~SmoothDiffuse() { }

	void configure() {
		BSDF::configure();
		/* Verify the input parameter and fix them if necessary */
		m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_reflectance->getValue(its);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0 
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);
			
		return m_reflectance->getValue(bRec.its)
			* (INV_PI * Frame::cosTheta(bRec.wo));
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0 
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0) 
			return Spectrum(0.0f);

		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		return m_reflectance->getValue(bRec.its);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		
		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = Frame::cosTheta(bRec.wo) * INV_PI;
		return m_reflectance->getValue(bRec.its) 
			* (INV_PI * Frame::cosTheta(bRec.wo));
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) 
				&& (name == "reflectance" || name == "diffuseReflectance")) {
			m_reflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_reflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_reflectance.get());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothDiffuse[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  reflectance = " << indent(m_reflectance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
};

// ================ Hardware shader implementation ================ 

class SmoothDiffuseShader : public Shader {
public:
	SmoothDiffuseShader(Renderer *renderer, const Texture *reflectance) 
		: Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_reflectance.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * 0.31831 * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<Shader> m_reflectanceShader;
};

Shader *SmoothDiffuse::createShader(Renderer *renderer) const { 
	return new SmoothDiffuseShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(SmoothDiffuseShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothDiffuse, "Smooth diffuse BRDF")
MTS_NAMESPACE_END
