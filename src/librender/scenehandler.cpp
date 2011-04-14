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

#include <mitsuba/core/platform.h>
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/scene.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

#define XMLLog(level, fmt, ...) Thread::getThread()->getLogger()->log(\
	level, NULL, __FILE__, __LINE__, "Near file offset %i: " fmt, \
	(int) m_parser->getSrcOffset(), ## __VA_ARGS__)

SceneHandler::SceneHandler(const SAXParser *parser,
	const ParameterMap &params, NamedObjectMap *namedObjects,
	bool isIncludedFile) : m_parser(parser), m_params(params),
		m_namedObjects(namedObjects), m_isIncludedFile(isIncludedFile) {
		m_pluginManager = PluginManager::getInstance();

	if (m_isIncludedFile) {
		SAssert(namedObjects != NULL);
	} else {
		SAssert(namedObjects == NULL);
		m_namedObjects = new NamedObjectMap();
	}

#if !defined(WIN32)
	setlocale(LC_NUMERIC, "C");
#endif
}

SceneHandler::~SceneHandler() {
	clear();
	if (!m_isIncludedFile)
		delete m_namedObjects;
}

void SceneHandler::clear() {
	if (!m_isIncludedFile) {
		for (NamedObjectMap::iterator it = m_namedObjects->begin();
				it != m_namedObjects->end(); ++it)
			if (it->second)
				it->second->decRef();
		m_namedObjects->clear();
	}
}

// -----------------------------------------------------------------------
//  Implementation of the SAX DocumentHandler interface
// -----------------------------------------------------------------------

void SceneHandler::startDocument() {
	clear();
}

void SceneHandler::endDocument() {
	SAssert(m_scene != NULL);
}

void SceneHandler::characters(const XMLCh* const name,
	const unsigned int length) {
}

Float SceneHandler::parseFloat(const std::string &name,
		const std::string &str, Float defVal) const {
	char *end_ptr = NULL;
	if (str == "") {
		if (defVal == -1)
			XMLLog(EError, "Missing floating point value (in <%s>)", name.c_str());
		return defVal;
	}
	Float result = (Float) std::strtod(str.c_str(), &end_ptr);
	if (*end_ptr != '\0')
		XMLLog(EError, "Invalid floating point value specified (in <%s>)", name.c_str());
	return result;
}

void SceneHandler::startElement(const XMLCh* const xmlName,
	AttributeList &xmlAttributes) {
	std::string name = transcode(xmlName);

	ParseContext context((name == "scene") ? NULL : &m_context.top());

	/* Convert attributes to ISO-8859-1 */
	for (unsigned int i=0; i<xmlAttributes.getLength(); i++) {
		std::string attrValue = transcode(xmlAttributes.getValue(i));
		if (attrValue.length() > 0 && attrValue.find('$') != attrValue.npos) {
			for (std::map<std::string, std::string>::const_iterator it = m_params.begin();
				it != m_params.end(); ++it) {
				std::string::size_type pos = 0;
				std::string searchString = "$" + it->first;
				while ((pos = attrValue.find(searchString, pos)) != std::string::npos) {
					attrValue.replace(pos, searchString.size(), it->second);
					++pos;
				}
			}
			if (attrValue.find('$') != attrValue.npos)
				XMLLog(EError, "The scene referenced an undefined parameter: \"%s\"", attrValue.c_str());
		}

		context.attributes[transcode(xmlAttributes.getName(i))] = attrValue;
	}
	if (name == "transform")
		m_transform = Transform();

	m_context.push(context);
}

void SceneHandler::endElement(const XMLCh* const xmlName) {
	std::string name = transcode(xmlName);
	ParseContext &context = m_context.top();
	std::string type = boost::to_lower_copy(context.attributes["type"]);
	context.properties.setPluginName(type);
	if (context.attributes.find("id") != context.attributes.end())
		context.properties.setID(context.attributes["id"]);

	ref<ConfigurableObject> object = NULL;

	/* Construct configurable objects */
	if (name == "scene") {
		object = m_scene = new Scene(context.properties);
	} else if (name == "shape") {
		object = static_cast<Shape *> (m_pluginManager->createObject(
			MTS_CLASS(Shape), context.properties));
	} else if (name == "sampler") {
		object = static_cast<Sampler *> (m_pluginManager->createObject(
			MTS_CLASS(Sampler), context.properties));
	} else if (name == "film") {
		object = static_cast<Film *> (m_pluginManager->createObject(
			MTS_CLASS(Film), context.properties));
	} else if (name == "integrator") {
		object = static_cast<Integrator *> (m_pluginManager->createObject(
			MTS_CLASS(Integrator), context.properties));
	} else if (name == "texture") {
		object = static_cast<Texture *> (m_pluginManager->createObject(
			MTS_CLASS(Texture), context.properties));
	} else if (name == "camera") {
		object = static_cast<Camera *> (m_pluginManager->createObject(
			MTS_CLASS(Camera), context.properties));
	} else if (name == "subsurface") {
		object = static_cast<Subsurface *> (m_pluginManager->createObject(
			MTS_CLASS(Subsurface), context.properties));
	} else if (name == "luminaire") {
		object = static_cast<Luminaire *> (m_pluginManager->createObject(
			MTS_CLASS(Luminaire), context.properties));
	} else if (name == "medium") {
		object = static_cast<Medium *> (m_pluginManager->createObject(
			MTS_CLASS(Medium), context.properties));
	} else if (name == "volume") {
		object = static_cast<VolumeDataSource *> (m_pluginManager->createObject(
			MTS_CLASS(VolumeDataSource), context.properties));
	} else if (name == "phase") {
		object = static_cast<PhaseFunction *> (m_pluginManager->createObject(
			MTS_CLASS(PhaseFunction), context.properties));
	} else if (name == "bsdf") {
		object = static_cast<BSDF *> (m_pluginManager->createObject(
			MTS_CLASS(BSDF), context.properties));
	} else if (name == "rfilter") {
		object = static_cast<ReconstructionFilter *> (m_pluginManager->createObject(
			MTS_CLASS(ReconstructionFilter), context.properties));
	} else if (name == "null") {
		object = NULL;
	} else if (name == "ref") {
		std::string id = context.attributes["id"];
		if (m_namedObjects->find(id) == m_namedObjects->end())
			XMLLog(EError, "Referenced object '%s' not found!", id.c_str());
		object = (*m_namedObjects)[id];
	/* Construct properties */
	} else if (name == "integer") {
		char *end_ptr = NULL;
#ifdef WIN32
		int64_t i = _strtoi64(context.attributes["value"].c_str(), &end_ptr, 10);
#else
		int64_t i = strtoll(context.attributes["value"].c_str(), &end_ptr, 10);
#endif
		if (*end_ptr != '\0')
			XMLLog(EError, "Invalid integer value specified (in <%s>)", 
				context.attributes["name"].c_str());
		context.parent->properties.setLong(context.attributes["name"], i);
	} else if (name == "float") {
		Float f = parseFloat(name, context.attributes["value"]);
		context.parent->properties.setFloat(context.attributes["name"], f);
	} else if (name == "boolean") {
		bool value = false;
		if (context.attributes["value"] == "true")
			value = true;
		context.parent->properties.setBoolean(context.attributes["name"],
			value);
	} else if (name == "string") {
		context.parent->properties.setString(context.attributes["name"],
			context.attributes["value"]);
	} else if (name == "translate") {
		Float x = parseFloat(name, context.attributes["x"], 0);
		Float y = parseFloat(name, context.attributes["y"], 0);
		Float z = parseFloat(name, context.attributes["z"], 0);
		m_transform = Transform::translate(Vector(x, y, z)) * m_transform;
	} else if (name == "rotate") {
		Float x = parseFloat(name, context.attributes["x"], 0);
		Float y = parseFloat(name, context.attributes["y"], 0);
		Float z = parseFloat(name, context.attributes["z"], 0);
		Float angle = parseFloat(name, context.attributes["angle"]);
		m_transform = Transform::rotate(Vector(x, y, z), angle) * m_transform;
	} else if (name == "lookAt") {
		Float ox = parseFloat(name, context.attributes["ox"]);
		Float oy = parseFloat(name, context.attributes["oy"]);
		Float oz = parseFloat(name, context.attributes["oz"]);
		Float tx = parseFloat(name, context.attributes["tx"]);
		Float ty = parseFloat(name, context.attributes["ty"]);
		Float tz = parseFloat(name, context.attributes["tz"]);
		Float ux = parseFloat(name, context.attributes["ux"], 0);
		Float uy = parseFloat(name, context.attributes["uy"], 0);
		Float uz = parseFloat(name, context.attributes["uz"], 0);
		Point o(ox, oy, oz), t(tx, ty, tz);
		Vector u(ux, uy, uz);

		if (u.lengthSquared() == 0) {
			/* If 'up' was not specified, use an arbitrary axis */
			Vector unused;
			coordinateSystem(normalize(t-o), u, unused);
		}

		m_transform =  Transform::lookAt(o, t, u) * m_transform;
	} else if (name == "scale") {
		bool hasXYZ = 
			context.attributes["x"] != "" ||
			context.attributes["y"] != "" ||
			context.attributes["z"] != "";
		bool hasValue = 
			context.attributes["value"] != "";
		Float x=0, y=0, z=0;

		if (hasXYZ && hasValue) {
			XMLLog(EError, "<scale>: provided both xyz and value arguments!");
		} else if (hasXYZ) {
			x = parseFloat(name, context.attributes["x"], 1);
			y = parseFloat(name, context.attributes["y"], 1);
			z = parseFloat(name, context.attributes["z"], 1);
		} else if (hasValue) {
			x = y = z = parseFloat(name, context.attributes["value"]);
		} else {
			XMLLog(EError, "<scale>: provided neither xyz nor value arguments!");
		}

		m_transform = Transform::scale(Vector(x, y, z)) * m_transform;
	} else if (name == "matrix") {
		std::vector<std::string> tokens = tokenize(
			context.attributes["value"], ", ");
		if (tokens.size() != 16)
			XMLLog(EError, "Invalid matrix specified");
		int index = 0;
		Matrix4x4 mtx;

		for (int i=0; i<4; ++i)
			for (int j=0; j<4; ++j)
				mtx.m[i][j] = parseFloat(name, tokens[index++]);

		m_transform = Transform(mtx) * m_transform;
	} else if (name == "point" || name == "vector") {
		Float x = parseFloat(name, context.attributes["x"]);
		Float y = parseFloat(name, context.attributes["y"]);
		Float z = parseFloat(name, context.attributes["z"]);

		context.parent->properties.setPoint(context.attributes["name"], Point(x, y, z));
	} else if (name == "rgb") {
		std::string valueStr = context.attributes["value"];
		std::vector<std::string> tokens = tokenize(valueStr, ", ");
		Float value[3];
		if (tokens.size() == 1 && tokens[0].length() == 7 && tokens[0][0] == '#') {
			char *end_ptr = NULL;
			/* Parse HTML-style hexadecimal colors */
			int encoded = strtol(tokens[0].c_str()+1, &end_ptr, 16);
			if (*end_ptr != '\0')
				XMLLog(EError, "Invalid rgb value specified (in <%s>)", context.attributes["name"].c_str());
			value[0] = ((encoded & 0xFF0000) >> 16) / 255.0f;
			value[1] = ((encoded & 0x00FF00) >> 8) / 255.0f;
			value[2] =  (encoded & 0x0000FF) / 255.0f;
		} else if (tokens.size() == 1) {
			value[0] = value[1] = value[2] = parseFloat(name, tokens[0]);
		} else if (tokens.size() == 3) {
			for (int i=0; i<3; i++) 
				value[i] = parseFloat(name, tokens[i]);
		} else {
			value[0] = value[1] = value[2] = 0; // avoid warning
			XMLLog(EError, "Invalid RGB value specified");
		}
		Spectrum specValue;
		specValue.fromLinearRGB(value[0], value[1], value[2]);
		context.parent->properties.setSpectrum(context.attributes["name"],
			specValue);
	} else if (name == "srgb") {
		std::string valueStr = context.attributes["value"];
		std::vector<std::string> tokens = tokenize(valueStr, ", ");
		Float value[3];
		if (tokens.size() == 1 && tokens[0].length() == 7 && tokens[0][0] == '#') {
			char *end_ptr = NULL;
			/* Parse HTML-style hexadecimal colors */
			int encoded = strtol(tokens[0].c_str()+1, &end_ptr, 16);
			if (*end_ptr != '\0')
				XMLLog(EError, "Invalid sRGB value specified (in <%s>)", context.attributes["name"].c_str());
			value[0] = ((encoded & 0xFF0000) >> 16) / 255.0f;
			value[1] = ((encoded & 0x00FF00) >> 8) / 255.0f;
			value[2] =  (encoded & 0x0000FF) / 255.0f;
		} else if (tokens.size() == 1) {
			value[0] = value[1] = value[2] = parseFloat(name, tokens[0]);
		} else if (tokens.size() == 3) {
			for (int i=0; i<3; i++) 
				value[i] = parseFloat(name, tokens[i]);
		} else {
			value[0] = value[1] = value[2] = 0; // avoid warning
			XMLLog(EError, "Invalid sRGB value specified");
		}
		Spectrum specValue;
		specValue.fromSRGB(value[0], value[1], value[2]);
		context.parent->properties.setSpectrum(context.attributes["name"],
			specValue);
	} else if (name == "blackbody") {
		Float temperature = parseFloat(name, context.attributes["temperature"]);
		BlackBodySpectrum bb(temperature);
		Spectrum discrete;
		discrete.fromSmoothSpectrum(&bb);
		context.parent->properties.setSpectrum(context.attributes["name"], discrete);
	} else if (name == "spectrum") {
		std::vector<std::string> tokens = tokenize(
			context.attributes["value"], ", ");
		Float value[SPECTRUM_SAMPLES];
		if (tokens.size() == 1) {
			value[0] = parseFloat(name, tokens[0]);
			context.parent->properties.setSpectrum(context.attributes["name"],
				Spectrum(value[0]));
		} else {
			if (tokens[0].find(':') != std::string::npos) {
				InterpolatedSpectrum interp(tokens.size());
				/* Wavelength -> Value mapping */
				for (size_t i=0; i<tokens.size(); i++) {
					std::vector<std::string> tokens2 = tokenize(tokens[i], ":");
					if (tokens2.size() != 2) 
						XMLLog(EError, "Invalid spectrum->value mapping specified");
					Float wavelength = parseFloat(name, tokens2[0]);
					Float value = parseFloat(name, tokens2[1]);
					interp.appendSample(wavelength, value);
				}
				Spectrum discrete;
				discrete.fromSmoothSpectrum(&interp);
				context.parent->properties.setSpectrum(context.attributes["name"],
					discrete);
			} else {
				if (tokens.size() != SPECTRUM_SAMPLES)
					XMLLog(EError, "Invalid spectrum value specified (incorrect length)");
				for (int i=0; i<SPECTRUM_SAMPLES; i++) 
					value[i] = parseFloat(name, tokens[i]);
				context.parent->properties.setSpectrum(context.attributes["name"],
					Spectrum(value));
			}
		}
	} else if (name == "transform") {
		context.parent->properties.setTransform(context.attributes["name"],
			m_transform);
		/* Do nothing */
	} else if (name == "include") {
		SAXParser* parser = new SAXParser();
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		fs::path schemaPath = resolver->resolveAbsolute("schema/scene.xsd");

		/* Check against the 'scene.xsd' XML Schema */
		parser->setDoSchema(true);
		parser->setValidationSchemaFullChecking(true);
		parser->setValidationScheme(SAXParser::Val_Always);
		parser->setExternalNoNamespaceSchemaLocation(schemaPath.file_string().c_str());

		/* Set the handler and start parsing */
		SceneHandler *handler = new SceneHandler(parser, m_params, m_namedObjects, true);
		parser->setDoNamespaces(true);
		parser->setDocumentHandler(handler);
		parser->setErrorHandler(handler);
		parser->setCalculateSrcOfs(true);
		fs::path path = resolver->resolve(context.attributes["filename"]);
		XMLLog(EInfo, "Parsing included file \"%s\" ..", path.filename().c_str());
		parser->parse(path.file_string().c_str());

		object = handler->getScene();
		delete parser;
		delete handler;
	} else {
		XMLLog(EError, "Unhandled tag \"%s\" encountered!", name.c_str());
	}

	if (object != NULL || name == "null") {
		std::string id = context.attributes["id"];
		std::string nodeName = context.attributes["name"];

		if (id != "" && name != "ref") {
			if (m_namedObjects->find(id) != m_namedObjects->end())
				XMLLog(EError, "Duplicate ID '%s' used in scene description!", id.c_str());
			(*m_namedObjects)[id] = object;
			if (object)
				object->incRef();
		}

		if (object) {
			/* If the object has a parent, add it to the parent's children list */
			if (context.parent != NULL) {
				object->incRef();
				context.parent->children.push_back(
					std::pair<std::string, ConfigurableObject *>(nodeName, object));
			}

			/* If the object has children, append them */
			for (std::vector<std::pair<std::string, ConfigurableObject *> >
					::iterator it = context.children.begin();
					it != context.children.end(); ++it) {
				if (it->second != NULL) {
					object->addChild(it->first, it->second);
					it->second->setParent(object);
					it->second->decRef();
				}
			}

			/* Don't configure a scene object if it is from an included file */
			if (name != "include" && (!m_isIncludedFile || !object->getClass()->derivesFrom(MTS_CLASS(Scene)))) 
				object->configure();
		}
	}

	/* Warn about unqueried properties */
	std::vector<std::string> unq = context.properties.getUnqueried();
	for (unsigned int i=0; i<unq.size(); ++i)
		XMLLog(EWarn, "Unqueried attribute \"%s\" in element \"%s\"", unq[i].c_str(), name.c_str());

	m_context.pop();
}

// -----------------------------------------------------------------------
//  Implementation of the SAX ErrorHandler interface
// -----------------------------------------------------------------------

void SceneHandler::warning(const SAXParseException& e) {
	SLog(EWarn, "Warning in file \"%s\" (line %i): %s",
		transcode(e.getSystemId()).c_str(), e.getLineNumber(),
		transcode(e.getMessage()).c_str());
}

void SceneHandler::error(const SAXParseException& e) {
	SLog(EError, "Error in file \"%s\" (line %i): %s",
		transcode(e.getSystemId()).c_str(), e.getLineNumber(),
		transcode(e.getMessage()).c_str());
}

void SceneHandler::fatalError(const SAXParseException& e) {
	SLog(EError, "Fatal error in file \"%s\" (line %i): %s",
		transcode(e.getSystemId()).c_str(), e.getLineNumber(),
		transcode(e.getMessage()).c_str());
}

MTS_NAMESPACE_END
