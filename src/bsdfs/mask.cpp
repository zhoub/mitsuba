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
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{mask}{Opacity mask}
 * \parameters{
 *     \parameter{opacity}{\Spectrum\Or\Texture}{
 *          Specifies the per-channel opacity (where $1=$ completely transparent)\default{0.5}. 
 *     }
 * }
 * \renderings{
 *     \rendering{Rendering without an opacity mask}
 *         {bsdf_mask_before.jpg}
 *     \rendering{Rendering \emph{with} an opacity mask (\lstref{mask-leaf})}
 *         {bsdf_mask_after.jpg}
 * }
 * This plugin applies an opacity mask to add nested BSDF instance. It interpolates
 * between perfectly transparent and completely opaque based on the \code{opacity}
 * parameter.
 *
 * The transparency is implemented as a forward-facing Diract delta distribution.
 * \vspace{5mm}
 *
 * \begin{xml}[caption=Material configuration for a transparent leaf,
 *    label=lst:mask-leaf]
 * <bsdf type="mask">
 *     <bsdf type="twosided">
 *         <bsdf type="diffuse">
 *             <texture name="reflectance" type="bitmap">
 *                 <string name="filename" value="leaf.jpg"/>
 *             </texture>
 *         </bsdf>
 *     </bsdf>
 *     <texture name="opacity" type="bitmap">
 *         <string name="filename" value="leaf_opacity.jpg"/>
 *         <float name="gamma" value="1"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */

class Mask : public BSDF {
public:
	Mask(const Properties &props) 
		: BSDF(props) {
		m_opacity = new ConstantSpectrumTexture(
		props.getSpectrum("opacity", Spectrum(0.5f)));
	}

	Mask(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_opacity = static_cast<Texture *>(manager->getInstance(stream));
		m_nestedBSDF = static_cast<BSDF *>(manager->getInstance(stream));
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_opacity.get());
		manager->serialize(stream, m_nestedBSDF.get());
	}

	void configure() {
		if (!m_nestedBSDF)
			Log(EError, "A child BSDF is required");
		m_components.clear();
		for (int i=0; i<m_nestedBSDF->getComponentCount(); ++i)
			m_components.push_back(m_nestedBSDF->getType(i));
		m_components.push_back(EDeltaTransmission | EFrontSide | EBackSide);
		m_usesRayDifferentials = m_nestedBSDF->usesRayDifferentials();
		m_opacity = ensureEnergyConservation(m_opacity, "opacity", 1.0f);
		BSDF::configure();
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		Spectrum opacity = m_opacity->getValue(bRec.its);

		if (measure == ESolidAngle)
			return m_nestedBSDF->eval(bRec, ESolidAngle) * opacity;
		else if (measure == EDiscrete && std::abs(1-dot(bRec.wi, -bRec.wo)) < Epsilon)
			return Spectrum(1.0f) - opacity;
		else
			return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		Float prob = m_opacity->getValue(bRec.its).getLuminance();

		if (measure == ESolidAngle)
			return m_nestedBSDF->pdf(bRec, ESolidAngle) * prob;
		else if (measure == EDiscrete && std::abs(1-dot(bRec.wi, -bRec.wo)) < Epsilon)
			return 1-prob;
		else
			return 0.0f;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &_sample) const {
		Point2 sample(_sample);
		Spectrum result(0.0f);

		Spectrum opacity = m_opacity->getValue(bRec.its);
		Float prob = opacity.getLuminance();

		bool sampleTransmission = bRec.typeMask & EDeltaTransmission
			&& (bRec.sampledComponent == -1 ||
				bRec.sampledComponent == getComponentCount()-1);
		bool sampleNested = bRec.sampledComponent == -1 || 
			bRec.sampledComponent < getComponentCount()-1;

		if (sampleTransmission && sampleNested) {
			if (sample.x <= prob) {
				sample.x /= prob;
				result = m_nestedBSDF->sample(bRec, pdf, sample);
				pdf *= prob;
			} else {
				bRec.wo = -bRec.wi;
				bRec.sampledComponent = getComponentCount()-1;
				bRec.sampledType = EDeltaTransmission;
				pdf = 1-prob;
				result = Spectrum(1.0f) - opacity;
			}
		} else if (sampleTransmission) {
			bRec.wo = -bRec.wi;
			bRec.sampledComponent = getComponentCount()-1;
			bRec.sampledType = EDeltaTransmission;
			pdf = 1;
			result = Spectrum(1.0f) - opacity;
		} else if (sampleNested) {
			result = m_nestedBSDF->sample(bRec, pdf, sample);
		}

		return result;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &_sample) const {
		Point2 sample(_sample);
		Spectrum opacity = m_opacity->getValue(bRec.its);
		Float prob = opacity.getLuminance();

		bool sampleTransmission = bRec.typeMask & EDeltaTransmission
			&& (bRec.sampledComponent == -1 ||
				bRec.sampledComponent == getComponentCount()-1);
		bool sampleNested = bRec.sampledComponent == -1 || 
			bRec.sampledComponent < getComponentCount()-1;

		if (sampleTransmission && sampleNested) {
			if (sample.x <= prob) {
				Float invProb = 1.0f / prob;
				sample.x *= invProb;
				return m_nestedBSDF->sample(bRec, sample) * invProb;
			} else {
				bRec.wo = -bRec.wi;
				bRec.sampledComponent = getComponentCount()-1;
				bRec.sampledType = EDeltaTransmission;
				return (Spectrum(1.0f) - opacity) / (1-prob);
			}
		} else if (sampleTransmission) {
			bRec.wo = -bRec.wi;
			bRec.sampledComponent = getComponentCount()-1;
			bRec.sampledType = EDeltaTransmission;
			return Spectrum(1.0f) - opacity;
		} else if (sampleNested) {
			return m_nestedBSDF->sample(bRec, sample);
		} else {
			return Spectrum(0.0f);
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "opacity") 
			m_opacity = static_cast<Texture *>(child);
		else if (child->getClass()->derivesFrom(MTS_CLASS(BSDF)))
			m_nestedBSDF = static_cast<BSDF *>(child);
		else
			BSDF::addChild(name, child);
	}

	Shader *createShader(Renderer *renderer) const;

	std::string toString() const {
		std::ostringstream oss;
		oss << "Mask[" << endl
			<< "  opacity = " << indent(m_opacity->toString()) << "," << endl
			<< "  nestedBSDF = " << indent(m_nestedBSDF.toString()) << endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
protected:
	ref<Texture> m_opacity;
	ref<BSDF> m_nestedBSDF;
};

// ================ Hardware shader implementation ================ 

/**
 * Somewhat lame GLSL version, which doesn't actually render
 * the object as transparent and only modulates the nested
 * BRDF reflectance by the opacity mask
 */
class MaskShader : public Shader {
public:
	MaskShader(Renderer *renderer, const Texture *opacity, const BSDF *bsdf) 
		: Shader(renderer, EBSDFShader), m_opacity(opacity), m_bsdf(bsdf) {
		m_opacityShader = renderer->registerShaderForResource(opacity);
		m_bsdfShader = renderer->registerShaderForResource(bsdf);
	}

	bool isComplete() const {
		return m_opacityShader.get() != NULL && m_bsdfShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_opacity);
		renderer->unregisterShaderForResource(m_bsdf);
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_opacityShader);
		deps.push_back(m_bsdfShader);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << depNames[0] << "(uv) * " << depNames[1] << "(uv, wi, wo);" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << depNames[0] << "(uv) * " << depNames[1] << "_diffuse(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_opacity;
	ref<Shader> m_opacityShader;
	ref<const BSDF> m_bsdf;
	ref<Shader> m_bsdfShader;
};

Shader *Mask::createShader(Renderer *renderer) const { 
	return new MaskShader(renderer, m_opacity.get(), m_nestedBSDF.get());
}

MTS_IMPLEMENT_CLASS(MaskShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Mask, false, BSDF)
MTS_EXPORT_PLUGIN(Mask, "Mask BSDF");
MTS_NAMESPACE_END
