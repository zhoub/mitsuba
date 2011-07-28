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

#include <mitsuba/core/properties.h>
#include <mitsuba/core/netobject.h>

MTS_NAMESPACE_BEGIN

bool Properties::hasProperty(const std::string &name) const {
	return m_elements.find(name) != m_elements.end();
}

std::vector<std::string> Properties::getUnqueried() const {
	std::map<std::string, Element>::const_iterator it = m_elements.begin();
	std::vector<std::string> result;

	for (; it != m_elements.end(); ++it) {
		if ((*it).second.queried == false) 
			result.push_back((*it).first);
	}

	return result;
}

Properties::EPropertyType Properties::getType(const std::string &name) const {
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if (it == m_elements.end())
		SLog(EError, "Property \"%s\" has not been specified!", name.c_str());
	return (*it).second.type;
}

void Properties::setBoolean(const std::string &name, bool value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EBoolean;
	m_elements[name].v_boolean = value;
	m_elements[name].queried = false;
}

bool Properties::getBoolean(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EBoolean)
		SLog(EError, "The property \"%s\" has the wrong type (expected <boolean>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_boolean;
}

bool Properties::getBoolean(const std::string &name, bool defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EBoolean)
		SLog(EError, "The property \"%s\" has the wrong type (expected <boolean>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_boolean;
}

void Properties::setInteger(const std::string &name, int value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EInteger;
	m_elements[name].v_long = value;
	m_elements[name].queried = false;
}

int Properties::getInteger(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EInteger)
		SLog(EError, "The property \"%s\" has the wrong type (expected <integer>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (int) (*it).second.v_long;
}

int Properties::getInteger(const std::string &name, int defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EInteger)
		SLog(EError, "The property \"%s\" has the wrong type (expected <integer>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (int) (*it).second.v_long;
}

void Properties::setLong(const std::string &name, int64_t value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EInteger;
	m_elements[name].v_long = value;
	m_elements[name].queried = false;
}

int64_t Properties::getLong(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EInteger)
		SLog(EError, "The property \"%s\" has the wrong type (expected <integer>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_long;
}

int64_t Properties::getLong(const std::string &name, int64_t defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EInteger)
		SLog(EError, "The property \"%s\" has the wrong type (expected <integer>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_long;
}

size_t Properties::getSize(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EInteger)
		SLog(EError, "The property \"%s\" has the wrong type (expected <integer>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	if ((*it).second.v_long < 0)
		SLog(EError, "Size property \"%s\": expected a nonnegative value!");
	(*it).second.queried = true;
	return (size_t) (*it).second.v_long;
}

size_t Properties::getSize(const std::string &name, size_t defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EInteger)
		SLog(EError, "The property \"%s\" has the wrong type (expected <integer>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	if ((*it).second.v_long < 0)
		SLog(EError, "Size property \"%s\": expected a nonnegative value!");
	(*it).second.queried = true;
	return (size_t) (*it).second.v_long;
}

void Properties::setData(const std::string &name, Data value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EData;
	m_elements[name].v_data = value;
	m_elements[name].queried = false;
}

Properties::Data Properties::getData(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EData)
		SLog(EError, "The property \"%s\" has the wrong type (expected <data>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_data;
}

Properties::Data Properties::getData(const std::string &name, Data defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EData)
		SLog(EError, "The property \"%s\" has the wrong type (expected <data>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_data;
}

void Properties::setFloat(const std::string &name, Float value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EFloat;
	m_elements[name].v_float = value;
	m_elements[name].queried = false;
}

Float Properties::getFloat(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EFloat)
		SLog(EError, "The property \"%s\" has the wrong type (expected <float>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_float;
}

Float Properties::getFloat(const std::string &name, Float defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EFloat)
		SLog(EError, "The property \"%s\" has the wrong type (expected <float>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_float;
}

void Properties::setTransform(const std::string &name, const Transform &value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = ETransform;
	m_elements[name].v_transform = value;
	m_elements[name].queried = false;
}

Transform Properties::getTransform(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != ETransform)
		SLog(EError, "The property \"%s\" has the wrong type (expected <transform>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_transform;
}

Transform Properties::getTransform(const std::string &name, const Transform &defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != ETransform)
		SLog(EError, "The property \"%s\" has the wrong type (expected <transform>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_transform;
}

void Properties::setSpectrum(const std::string &name, const Spectrum &value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = ESpectrum;
	m_elements[name].v_spectrum = value;
	m_elements[name].queried = false;
}

Spectrum Properties::getSpectrum(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != ESpectrum)
		SLog(EError, "The property \"%s\" has the wrong type (expected <spectrum>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_spectrum;
}

Spectrum Properties::getSpectrum(const std::string &name, const Spectrum &defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != ESpectrum)
		SLog(EError, "The property \"%s\" has the wrong type (expected <spectrum>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_spectrum;
}

void Properties::setString(const std::string &name, const std::string &value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EString;
	m_elements[name].v_string = value;
	m_elements[name].queried = false;
}

std::string Properties::getString(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EString)
		SLog(EError, "The property \"%s\" has the wrong type (expected <string>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_string;
}

std::string Properties::getString(const std::string &name, const std::string &defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EString)
		SLog(EError, "The property \"%s\" has the wrong type (expected <string>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_string;
}

void Properties::setPoint(const std::string &name, const Point &value, bool warnDuplicates) {
	if (hasProperty(name) && warnDuplicates)
		SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str());
	m_elements[name].type = EPoint;
	m_elements[name].v_point = value;
	m_elements[name].queried = false;
}

Vector Properties::getVector(const std::string &name) const {
	return Vector(getPoint(name));
}

Vector Properties::getVector(const std::string &name, const Vector &defVal) const {
	return Vector(getPoint(name, Point(defVal)));
}

Point Properties::getPoint(const std::string &name) const {
	if (!hasProperty(name))
		SLog(EError, "Property \"%s\" missing", name.c_str());
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EPoint)
		SLog(EError, "The property \"%s\" has the wrong type (expected <point>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_point;
}

Point Properties::getPoint(const std::string &name, const Point &defVal) const {
	if (!hasProperty(name))
		return defVal;
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if ((*it).second.type != EPoint)
		SLog(EError, "The property \"%s\" has the wrong type (expected <point>). The detailed "
				"listing was:\n%s", name.c_str(), toString().c_str());
	(*it).second.queried = true;
	return (*it).second.v_point;
}

std::string Properties::toString() const {
	std::map<std::string, Element>::const_iterator it = m_elements.begin();
	std::ostringstream oss;

	oss << "Properties[" << endl
		<< "  pluginName = \"" << m_pluginName << "\"," << endl
		<< "  elements = {" << endl;
	for (; it != m_elements.end(); ++it) {
		oss << "    \"" << (*it).first << "\" -> ";
		const Element &el = (*it).second;
		switch (el.type) {
			case EBoolean:
				oss << (el.v_boolean ? "true" : "false");
				break;
			case EInteger:
				oss << el.v_long;
				break;
			case EFloat:
				oss << el.v_float;
				break;
			case EPoint:
				oss << el.v_point.toString();
				break;
			case ETransform:
				oss << indent(el.v_transform.toString(), 2);
				break;
			case ESpectrum:
				oss << indent(el.v_spectrum.toString(), 2);
				break;
			case EString:
				oss << "\"" << el.v_string << "\"";
				break;
			case EData:
				oss << el.v_data.ptr << " (size=" << el.v_data.size << ")";
				break;
			default:
				SLog(EError, "Encountered an unknown property type!");
		}
		oss << "," << endl;
	}
	oss << "  }" << endl
		<< "]" << endl;
	return oss.str();
}

void Properties::markQueried(const std::string &name) const {
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if (it == m_elements.end())
		SLog(EError, "Could not find parameter \"%s\"!", name.c_str());
	it->second.queried = true;
}

bool Properties::wasQueried(const std::string &name) const {
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if (it == m_elements.end())
		SLog(EError, "Could not find parameter \"%s\"!", name.c_str());
	return it->second.queried;
}

ConfigurableObject::ConfigurableObject(Stream *stream, InstanceManager *manager) 
 : SerializableObject(stream, manager), m_configured(true) {
	m_parent = static_cast<ConfigurableObject *>(manager->getInstance(stream));
}

void ConfigurableObject::setParent(ConfigurableObject *parent) {
	m_parent = parent;
}

void ConfigurableObject::configure() {
}

void ConfigurableObject::serialize(Stream *stream, InstanceManager *manager) const {
	if (!getClass()->isSerializable())
		Log(EError, "Error: trying to serialize an instance of type '%s', which does "
			"not have full serialization support!", getClass()->getName().c_str());
	manager->serialize(stream, m_parent);
}

void ConfigurableObject::addChild(const std::string &name, ConfigurableObject *child) {
	SLog(EError, "ConfigurableObject::addChild(\"%s\", %s) not implemented in \"%s\"", 
		name.c_str(), child->toString().c_str(), toString().c_str());
}

void NetworkedObject::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
}
	
void NetworkedObject::bindUsedResources(ParallelProcess *proc) const {
}

void NetworkedObject::wakeup(std::map<std::string, SerializableObject *> &params) {
}

MTS_IMPLEMENT_CLASS(ConfigurableObject, true, SerializableObject)
MTS_IMPLEMENT_CLASS(NetworkedObject, true, ConfigurableObject)
MTS_NAMESPACE_END
