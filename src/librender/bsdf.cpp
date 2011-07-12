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
#include <mitsuba/core/frame.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

BSDF::BSDF(const Properties &props)
 : ConfigurableObject(props), m_name(props.getID()) {
	/* By default, verify whether energy conservation holds
	   for the user-specified parameter values. This step
	   is completely up to the particular BSDF implementations */
	m_ensureEnergyConservation = props.getBoolean(
		"ensureEnergyConservation", true);
}

BSDF::BSDF(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_name = stream->readString();
	m_ensureEnergyConservation = stream->readBool();
}

BSDF::~BSDF() { }

void BSDF::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	stream->writeString(m_name);
	stream->writeBool(m_ensureEnergyConservation);
}

void BSDF::setParent(ConfigurableObject *parent) {
	/* BSDF's don't need to reference their parent -> do nothing */
}

void BSDF::addChild(const std::string &name, ConfigurableObject *obj) {
	ConfigurableObject::addChild(name, obj);
}

void BSDF::configure() {
	m_combinedType = 0;
	for (size_t i=0; i<m_components.size(); ++i)
		m_combinedType |= m_components[i];
}

Spectrum BSDF::getDiffuseReflectance(const Intersection &its) const {
	BSDFQueryRecord bRec(its);
	bRec.wi = bRec.wo = Vector(0, 0, 1);
	bRec.typeMask = EDiffuseReflection;
	return eval(bRec) * M_PI;
}

Texture *BSDF::ensureEnergyConservation(Texture *texture, 
		const std::string &paramName, Float max) const {
	if (!m_ensureEnergyConservation)
		return texture;

	Float actualMax = texture->getMaximum().max();
	if (actualMax > max) {
		std::ostringstream oss;
		Float scale = 0.99f * (max / actualMax);
		oss << "The BSDF" << endl << toString() << endl
			<< "violates energy conservation! The parameter \"" << paramName << "\" " 
			<< "has a component-wise maximum of "<< actualMax << " (which is > " << max << "!) "
			<< "and will therefore be scaled by " << scale << " to prevent "
			<< "issues. Specify the parameter ensureEnergyConservation=false "
			<< "to the BSDF to prevent this from happening.";
		Log(EWarn, "%s", oss.str().c_str());
		Properties props("scale");
		props.setFloat("value", scale);
		Texture *scaleTexture = static_cast<Texture *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Texture), props));
		scaleTexture->addChild("", texture);
		scaleTexture->configure();
		return scaleTexture;
	}
	return texture;
}

std::pair<Texture *, Texture *> BSDF::ensureEnergyConservation(
		Texture *tex1, Texture *tex2, const std::string &paramName1,
		const std::string &paramName2, Float max) const {
	if (!m_ensureEnergyConservation)
		return std::make_pair(tex1, tex2);
	Float actualMax = (tex1->getMaximum() + tex2->getMaximum()).max();
	if (actualMax > max) {
		std::ostringstream oss;
		Float scale = 0.99f * (max / actualMax);
		oss << "The BSDF" << endl << toString() << endl
			<< "violates energy conservation! The parameters \"" << paramName1 << "\" " 
			<< "and \"" << paramName2 << "\" sum to a component-wise maximum of "
			<< actualMax << " (which is > " << max << "!) and will therefore be "
			<< "scaled by " << scale << " to prevent issues. Specify the parameter "
			<< "ensureEnergyConservation=false to the BSDF to prevent this from "
			<< "happening.";
		Log(EWarn, "%s", oss.str().c_str());
		Properties props("scale");
		props.setFloat("value", scale);
		Texture *scaleTexture1 = static_cast<Texture *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Texture), props));
		Texture *scaleTexture2 = static_cast<Texture *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Texture), props));
		scaleTexture1->addChild("", tex1);
		scaleTexture1->configure();
		scaleTexture2->addChild("", tex2);
		scaleTexture2->configure();
		return std::make_pair(scaleTexture1, scaleTexture2);
	}

	return std::make_pair(tex1, tex2);
}

static std::string typeMaskToString(unsigned int typeMask) {
	std::ostringstream oss;
	oss << "{ ";
	#define isset(mask) (typeMask & mask) == mask
	if (isset(BSDF::EAll)) { oss << "all "; typeMask &= ~BSDF::EAll; }
	if (isset(BSDF::ESmooth)) { oss << "smooth "; typeMask &= ~BSDF::ESmooth; }
	if (isset(BSDF::EDiffuse)) { oss << "diffuse "; typeMask &= ~BSDF::EDiffuse; }
	if (isset(BSDF::EGlossy)) { oss << "glossy "; typeMask &= ~BSDF::EGlossy; }
	if (isset(BSDF::EDelta)) { oss << "delta"; typeMask &= ~BSDF::EDelta; }
	if (isset(BSDF::EDelta1D)) { oss << "delta1D "; typeMask &= ~BSDF::EDelta1D; }
	if (isset(BSDF::EDiffuseReflection)) { oss << "diffuseReflection "; typeMask &= ~BSDF::EDiffuseReflection; }
	if (isset(BSDF::EDiffuseTransmission)) { oss << "diffuseTransmission "; typeMask &= ~BSDF::EDiffuseTransmission; }
	if (isset(BSDF::EGlossyReflection)) { oss << "glossyReflection "; typeMask &= ~BSDF::EGlossyReflection; }
	if (isset(BSDF::EGlossyTransmission)) { oss << "glossyTransmission "; typeMask &= ~BSDF::EGlossyTransmission; }
	if (isset(BSDF::EDeltaReflection)) { oss << "deltaReflection "; typeMask &= ~BSDF::EDeltaReflection; }
	if (isset(BSDF::EDeltaTransmission)) { oss << "deltaTransmission "; typeMask &= ~BSDF::EDeltaTransmission; }
	if (isset(BSDF::EDelta1DReflection)) { oss << "delta1DReflection "; typeMask &= ~BSDF::EDelta1DReflection; }
	if (isset(BSDF::EDelta1DTransmission)) { oss << "delta1DTransmission "; typeMask &= ~BSDF::EDelta1DTransmission; }
	if (isset(BSDF::EAnisotropic)) { oss << "anisotropic "; typeMask &= ~BSDF::EAnisotropic; }
	if (isset(BSDF::EFrontSide)) { oss << "frontSide "; typeMask &= ~BSDF::EFrontSide; }
	if (isset(BSDF::EBackSide)) { oss << "backSide "; typeMask &= ~BSDF::EBackSide; }
	if (isset(BSDF::ECanUseSampler)) { oss << "canUseSampler "; typeMask &= ~BSDF::ECanUseSampler; }
	#undef isset
	oss << "}";
	return oss.str();
}

std::string BSDFQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "BSDFQueryRecord[" << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  wo = " << wo.toString() << "," << std::endl
		<< "  quantity = " << quantity << "," << std::endl
		<< "  typeMask = " << typeMaskToString(typeMask) << "," << std::endl
		<< "  sampledType = " << typeMaskToString(sampledType) << "," << std::endl
		<< "  component = " << component << "," << std::endl
		<< "  sampledComponent = " << sampledComponent << "," << std::endl
		<< "]";
	return oss.str();
}
	
MTS_IMPLEMENT_CLASS(BSDF, true, ConfigurableObject)
MTS_NAMESPACE_END
