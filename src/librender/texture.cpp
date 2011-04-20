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

#include <mitsuba/render/consttexture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

Texture::Texture(const Properties &props)
 : ConfigurableObject(props) {
}

Texture::Texture(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
}
	
Vector3i Texture::getResolution() const {
	return Vector3i(0);
}

Texture::~Texture() {
}

void Texture::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
}

Texture2D::Texture2D(const Properties &props) : Texture(props) {
	if (props.getString("coordinates", "uv") == "uv") {
		m_uvOffset = Point2(
			props.getFloat("uoffset", 0.0f),
			props.getFloat("voffset", 0.0f)
		);
		m_uvScale = Vector2(
			props.getFloat("uscale", 1.0f),
			props.getFloat("vscale", 1.0f)
		);
	} else {
		Log(EError, "Only UV coordinates are supported at the moment!");
	}
}

Texture2D::Texture2D(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_uvOffset = Point2(stream);
	m_uvScale = Vector2(stream);
}

Texture2D::~Texture2D() {
}

void Texture2D::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);
	m_uvOffset.serialize(stream);
	m_uvScale.serialize(stream);
}

Spectrum Texture2D::getValue(const Intersection &its) const {
	Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;
	if (its.hasUVPartials) {
		return getValue(uv, 
			its.dudx * m_uvScale.x, its.dudy * m_uvScale.x,
			its.dvdx * m_uvScale.y, its.dvdy * m_uvScale.y);
	} else {
		return getValue(uv);
	}
}

ConstantTexture::ConstantTexture(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_value = Spectrum(stream);
}

void ConstantTexture::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);

	m_value.serialize(stream);
}

class ConstantTextureShader : public Shader {
public:
	ConstantTextureShader(Renderer *renderer, const Spectrum &value) 
		: Shader(renderer, ETextureShader), m_value(value) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_value;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return " << evalName << "_value;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_value"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_value);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_value;
};

Shader *ConstantTexture::createShader(Renderer *renderer) const { 
	return new ConstantTextureShader(renderer, m_value);
}

MTS_IMPLEMENT_CLASS(Texture, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(Texture2D, true, Texture)
MTS_IMPLEMENT_CLASS_S(ConstantTexture, false, Texture)
MTS_IMPLEMENT_CLASS(ConstantTextureShader, false, Shader)
MTS_NAMESPACE_END
