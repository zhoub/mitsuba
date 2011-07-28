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
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/basicshader.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{hk}{Hanrahan-Krueger BSDF}
 * \icon{bsdf_hk}
 *
 * \parameters{
 *     \parameter{material}{\String}{Name of a material preset, see 
 *           \tblref{medium-coefficients}.\!\default{\texttt{skin1}}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{Scattering coefficient of the 
 *      layer. \default{based on \code{material}}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{Absorption coefficient of the 
 *      layer. \default{based on \code{material}}}
 *     \parameter{thickness}{\Float}{Denotes the thickness of the layer.
 *      Should be specified in inverse units of \code{sigmaA} and \code{sigmaS})\default{1}}
 *     \parameter{\Unnamed}{\Phase}{A nested phase function instance that represents 
 *      the type of scattering interactions occurring within the layer}
 * }
 *
 * \renderings{
 *     \rendering{$\sigma_s=2$, $\sigma_a=1$, thickness$=0.1$}{bsdf_hk_1}
 *     \rendering{\code{ketchup} material preset}{bsdf_hk_2}
 * }
 *
 * This plugin provides an implementation of the Hanrahan-Krueger BSDF
 * \cite{Hanrahan1993Reflection} for simulating single scattering in thin 
 * index-matched layers filled with random scattering media. 
 * Apart from single scattering, the implementation also accounts for attenuated
 * light that passes through the medium without undergoing any scattering events.
 *
 * This BSDF requires a phase function to model scattering interactions within the
 * random medium. When no phase function is explicitly specified, it uses an 
 * isotropic one ($g=0$) by default. A sample usage for instantiating the
 * plugin is given below: 
 * \begin{xml}
 * <bsdf type="hk">
 *     <spectrum name="sigmaS" value="2"/>
 *     <spectrum name="sigmaA" value="0.05"/>
 *     <float name="thickness" value="1"/>
 *
 *     <phase type="hg">
 *         <float name="g" value="0.8"/>
 *     </phase>
 * </bsdf>
 * \end{xml}
 *
 * When used in conjuction with the \pluginref{coating} plugin, it is possible 
 * to model refraction and reflection at the layer boundaries when the indices 
 * of refraction are mismatched. The combination of these two plugins
 * reproduces the full model proposed in \cite{Hanrahan1993Reflection}.
 *
 * \begin{xml}
 * <bsdf type="coating">
 *     <float name="extIOR" value="1.0"/>
 *     <float name="intIOR" value="1.5"/>
 *
 *     <bsdf type="hk">
 *         <spectrum name="sigmaS" value="2"/>
 *         <spectrum name="sigmaA" value="0.05"/>
 *         <float name="thickness" value="1"/>
 *
 *         <phase type="hg">
 *             <float name="g" value="0.8"/>
 *         </phase>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 *
 * Note that when \texttt{sigmaS} = \texttt{sigmaA}$\ = 0$, or when \texttt{thickness=0},
 * any geometry associated with this scattering model will be invisible.
 *
 * The implementation in Mitsuba is based on code by Tom Kazimiers and Marios Papas.
*/

class HanrahanKrueger : public BSDF {
public:
	HanrahanKrueger(const Properties &props) : BSDF(props) {
		Spectrum sigmaS, sigmaA;
		lookupMaterial(props, sigmaS, sigmaA);

		/* Scattering coefficient of the layer */
		m_sigmaS = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaS", sigmaS));

		/* Absorption coefficient of the layer */
		m_sigmaA = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaA", sigmaA));

		/* Slab thickness in inverse units of sigmaS and sigmaA */
		m_thickness = props.getFloat("thickness", 1); 

	}

	HanrahanKrueger(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_phase = static_cast<PhaseFunction *>(manager->getInstance(stream));
		m_sigmaS = static_cast<Texture *>(manager->getInstance(stream));
		m_sigmaA = static_cast<Texture *>(manager->getInstance(stream));
		m_thickness = stream->readFloat();
		configure();
	}

	void configure() {
		if (m_phase == NULL)
			m_phase = static_cast<PhaseFunction *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(PhaseFunction), Properties("isotropic")));

		m_components.clear();
		m_components.push_back(EGlossyReflection   | EFrontSide | EBackSide | ECanUseSampler);
		m_components.push_back(EGlossyTransmission | EFrontSide | EBackSide | ECanUseSampler);
		m_components.push_back(EDeltaTransmission  | EFrontSide | EBackSide | ECanUseSampler);

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Spectrum sigmaA = m_sigmaA->getValue(its),
				 sigmaS = m_sigmaS->getValue(its),
				 sigmaT = sigmaA + sigmaS,
				 albedo;
		for (int i = 0; i < SPECTRUM_SAMPLES; i++)
			albedo[i] = sigmaT[i] > 0 ? (sigmaS[i]/sigmaT[i]) : (Float) 0;
		return albedo; /* Very approximate .. */
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		Spectrum sigmaA = m_sigmaA->getValue(bRec.its),
				 sigmaS = m_sigmaS->getValue(bRec.its),
				 sigmaT = sigmaA + sigmaS,
				 tauD = sigmaT * m_thickness,
				 result(0.0f);

		if (measure == EDiscrete) {
			/* Figure out if the specular transmission is specifically requested */
			bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 2);

			/* Return the attenuated light if requested */
			if (hasSpecularTransmission &&
				std::abs(1+dot(bRec.wi, bRec.wo)) < Epsilon)
				result = (-tauD/std::abs(Frame::cosTheta(bRec.wi))).exp();
		} else if (measure == ESolidAngle) {
			/* Sample single scattering events */
			bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
			bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		
			Spectrum albedo;
			for (int i = 0; i < SPECTRUM_SAMPLES; i++)
				albedo[i] = sigmaT[i] > 0 ? (sigmaS[i]/sigmaT[i]) : (Float) 0;


			const Float cosThetaI = Frame::cosTheta(bRec.wi),
				        cosThetaO = Frame::cosTheta(bRec.wo),
				        dp = cosThetaI*cosThetaO;
		
			bool reflection = dp > 0, transmission = dp < 0;

			/* ==================================================================== */
			/*                        Reflection component                          */
			/* ==================================================================== */

			if (hasGlossyReflection && reflection) {
				MediumSamplingRecord dummy;
				PhaseFunctionQueryRecord pRec(dummy,bRec.wi,bRec.wo); 
				const Float phaseVal = m_phase->eval(pRec);

				result = albedo * (phaseVal*cosThetaI/(cosThetaI+cosThetaO)) *
					(Spectrum(1.0f)-((-tauD/std::abs(cosThetaI))+(-tauD/std::abs(cosThetaO))).exp());
			}

			/* ==================================================================== */
			/*                       Transmission component                         */
			/* ==================================================================== */

			if (hasGlossyTransmission && transmission) {
				MediumSamplingRecord dummy;
				PhaseFunctionQueryRecord pRec(dummy,bRec.wi,bRec.wo);
				const Float phaseVal = m_phase->eval(pRec);

				/* Hanrahan etal 93 Single Scattering transmission term */
				if (std::abs(cosThetaI + cosThetaO) < Epsilon) {
					/* avoid division by zero */
					result += albedo * phaseVal*tauD/std::abs(cosThetaO) * 
								((-tauD/std::abs(cosThetaO)).exp());
				} else {
					/* Guaranteed to be positive even if |cosThetaO| > |cosThetaI| */
					result += albedo * phaseVal*std::abs(cosThetaI)/(std::abs(cosThetaI)-std::abs(cosThetaO)) * 
						((-tauD/std::abs(cosThetaI)).exp() - (-tauD/std::abs(cosThetaO)).exp());
				}
			}
		}


		return result;
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool hasSingleScattering = (bRec.typeMask & EGlossy)
			&& (bRec.component == -1 || bRec.component == 0 || bRec.component == 1);
		bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
			&& (bRec.component == -1 || bRec.component == 2);

		const Spectrum sigmaA = m_sigmaA->getValue(bRec.its),
				 sigmaS = m_sigmaS->getValue(bRec.its),
				 sigmaT = sigmaA + sigmaS,
				 tauD = sigmaT * m_thickness;

		Float probSpecularTransmission = (-tauD/std::abs(Frame::cosTheta(bRec.wi))).exp().average();

		if (measure == EDiscrete) {
			bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 2);
			/* Return the attenuated light if requested */
			if (hasSpecularTransmission &&
				std::abs(1+dot(bRec.wi, bRec.wo)) < Epsilon)
				return hasSingleScattering ? probSpecularTransmission : 1.0f;
		} else if (hasSingleScattering && measure == ESolidAngle) {
			bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
			bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
			bool reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0;

			if ((!hasGlossyReflection && reflection) ||
				(!hasGlossyTransmission && !reflection))
				return 0.0f;

			/* Sampled according to the phase function lobe(s) */
			MediumSamplingRecord dummy;
			PhaseFunctionQueryRecord pRec(dummy, bRec.wi, bRec.wo);
			Float pdf = m_phase->pdf(pRec);
			if (hasSpecularTransmission)
				pdf *= 1-probSpecularTransmission;
			return pdf;
		}
		return 0.0f;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, Float &_pdf, const Point2 &_sample) const {
		AssertEx(bRec.sampler != NULL, "The BSDFQueryRecord needs to have a sampler!");

		bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
			&& (bRec.component == -1 || bRec.component == 2);
		bool hasSingleScattering = (bRec.typeMask & EGlossy)
			&& (bRec.component == -1 || bRec.component == 0 || bRec.component == 1);

		const Spectrum sigmaA = m_sigmaA->getValue(bRec.its),
				 sigmaS = m_sigmaS->getValue(bRec.its),
				 sigmaT = sigmaA + sigmaS,
				 tauD = sigmaT * m_thickness;

		/* Probability for a specular transmission is approximated by the average (per wavelength) 
		 * probability of a photon exiting without a scattering event or an absorption event */
		Float probSpecularTransmission = (-tauD/std::abs(Frame::cosTheta(bRec.wi))).exp().average();

		bool choseSpecularTransmission = hasSpecularTransmission;

		Point2 sample(_sample);
		if (hasSpecularTransmission && hasSingleScattering) {
			if (sample.x > probSpecularTransmission) {
				sample.x = (sample.x - probSpecularTransmission) / (1 - probSpecularTransmission);
				choseSpecularTransmission = false;
			}
		}

		if (choseSpecularTransmission) {
			/* The specular transmission component was sampled */
			bRec.sampledComponent = 2;
			bRec.sampledType = EDeltaTransmission;

			bRec.wo = -bRec.wi;

			_pdf = hasSingleScattering ? probSpecularTransmission : 1.0f;
			return eval(bRec, EDiscrete);
		} else {
			/* The glossy transmission/scattering component should be sampled */
			bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
			bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

			/* Sample According to the phase function lobes */
			PhaseFunctionQueryRecord pRec(MediumSamplingRecord(), bRec.wi, bRec.wo);
			m_phase->sample(pRec, _pdf, bRec.sampler);

			/* Store the sampled direction */
			bRec.wo = pRec.wo;
			
			bool reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0;
			if ((!hasGlossyReflection && reflection) ||
				(!hasGlossyTransmission && !reflection))
				return Spectrum(0.0f);

			/* Notify that the scattering component was sampled */
			bRec.sampledComponent = reflection ? 0 : 1;
			bRec.sampledType = EGlossy;

			_pdf *= (hasSpecularTransmission ? (1 - probSpecularTransmission) : 1.0f);

			/* Guard against numerical imprecisions */
			if (_pdf == 0) 
				return Spectrum(0.0f);
			else
				return eval(bRec, ESolidAngle);

		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf = 0;
		Spectrum result = HanrahanKrueger::sample(bRec, pdf, sample);

		if (result.isZero())
			return Spectrum(0.0f);
		else
			return result / pdf;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_phase.get());
		manager->serialize(stream, m_sigmaS.get());
		manager->serialize(stream, m_sigmaA.get());
		stream->writeFloat(m_thickness);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
			Assert(m_phase == NULL);
			m_phase = static_cast<PhaseFunction *>(child);
		} else {
			Log(EError, "Invalid child node! (\"%s\")",
				cClass->getName().c_str());
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HanrahanKrueger[" << endl
   			<< "  sigmaS = " << indent(m_sigmaS->toString()) << "," << endl
   			<< "  sigmaA = " << indent(m_sigmaA->toString()) << "," << endl
   			<< "  phase = " << indent(m_phase->toString()) << "," << endl
   			<< "  thickness = " << m_thickness << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<PhaseFunction> m_phase;
	ref<Texture> m_sigmaS;
	ref<Texture> m_sigmaA;
	Float m_thickness;
};

MTS_IMPLEMENT_CLASS_S(HanrahanKrueger, false, BSDF)
MTS_EXPORT_PLUGIN(HanrahanKrueger, "Hanrahan-Krueger BSDF");
MTS_NAMESPACE_END
