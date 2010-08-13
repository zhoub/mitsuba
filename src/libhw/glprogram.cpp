#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/glprogram.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

extern GLEWContext *glewGetContext();

GLProgram::GLProgram(const std::string &name) 
 : GPUProgram(name) {
	m_id[0] = m_id[1] = m_program = 0;
}

void GLProgram::init() {
	Assert(m_id[0] == 0 && m_id[1] == 0 && m_program == 0);

	Log(EDebug, "Uploading a GPU program : %s", toString().c_str());

	m_id[EVertexProgram] = createShader(GL_VERTEX_SHADER_ARB, 
		m_source[EVertexProgram]);
	m_id[EFragmentProgram] = createShader(GL_FRAGMENT_SHADER_ARB,
		m_source[EFragmentProgram]);
	m_id[EGeometryProgram] = createShader(GL_GEOMETRY_SHADER_ARB,
		m_source[EGeometryProgram]);

	m_program = glCreateProgramObjectARB();
	if (m_id[EVertexProgram] != 0)
		glAttachObjectARB(m_program, m_id[EVertexProgram]);
	if (m_id[EFragmentProgram] != 0)
		glAttachObjectARB(m_program, m_id[EFragmentProgram]);
	if (m_id[EGeometryProgram] != 0)
		glAttachObjectARB(m_program, m_id[EGeometryProgram]);

	if (m_id[EGeometryProgram] != 0) {
		glProgramParameteriEXT(m_program, GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES); 
		glProgramParameteriEXT(m_program, GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
		glProgramParameteriEXT(m_program, GL_GEOMETRY_VERTICES_OUT_EXT, m_maxVertices);
	}

	glLinkProgramARB(m_program);
	//glValidateProgramARB(m_program);

	std::string infoLog = getInfoLogProgram();

	GLint result;
	glGetObjectParameterivARB(m_program, GL_LINK_STATUS, &result);
	if (result == GL_FALSE) {
		cleanup();
		if (infoLog != "")
			Log(EWarn, "%s", infoLog.c_str());
		Log(EError, "Error linking a GPU program!");
	} else if (infoLog != "") {
		if (infoLog.find("warning") != std::string::npos)
			Log(EWarn, infoLog.c_str());
		else
			Log(EDebug, infoLog.c_str());
	}
}

int GLProgram::createShader(int type, const std::string &source) {
	if (source == "")
		return 0;

	int id = glCreateShaderObjectARB(type);

	const char *string = source.c_str();
	GLint stringLength = (GLint) source.length();
	glShaderSourceARB(id, 1, &string, &stringLength);
	glCompileShaderARB(id);

	std::string infoLog = getInfoLogShader(id);

	GLint result;
	glGetObjectParameterivARB(id, GL_COMPILE_STATUS, &result);
	if (result == GL_FALSE) {
		cleanup();
		if (infoLog != "")
			Log(EError, "Error compiling a shader: %s", infoLog.c_str());
		else
			Log(EError, "Unknown error encountered while compiling a shader!");
	} else if (infoLog != "") {
		if (infoLog.find("warning") != std::string::npos)
			Log(EWarn, infoLog.c_str());
		else
			Log(EDebug, infoLog.c_str());
	}
	return id;
}

std::string GLProgram::getInfoLogShader(int id) {
	GLint stringLength;
	std::string result;

	glGetObjectParameterivARB(id, GL_INFO_LOG_LENGTH, &stringLength);
	if (stringLength > 0) {
		char *tmp = new char[stringLength + 1];
		glGetInfoLogARB(id, stringLength, &stringLength, tmp);
		result = tmp;
		delete[] tmp;
	}
	return result;
}

std::string GLProgram::getInfoLogProgram() {
	GLint stringLength;
	std::string result;

	glGetObjectParameterivARB(m_program, GL_INFO_LOG_LENGTH, &stringLength);
	if (stringLength > 0) {
		char *tmp = new char[stringLength + 1];
		glGetInfoLogARB(m_program, stringLength, &stringLength, tmp);
		result = tmp;
		delete[] tmp;
	}
	return result;
}

int GLProgram::getParameterID(const std::string &name, bool failIfMissing) const {
	int id = glGetUniformLocation(m_program, name.c_str());
	if (id == -1 && failIfMissing)
		Log(EError, "Unable to find the parameter named \"%s\"", name.c_str());
	return id;
}

void GLProgram::setParameter(int id, bool value) {
	if (id == -1)
		return;
	glUniform1i(id, value ? 1 : 0);
}

void GLProgram::setParameter(int id, Float value) {
	if (id == -1)
		return;
	glUniform1f(id, (GLfloat) value);
}

void GLProgram::setParameter(int id, const Vector &value) {
	if (id == -1)
		return;
	glUniform3f(id, (GLfloat) value.x, (GLfloat) value.y, (GLfloat) value.z);
}

void GLProgram::setParameter(int id, const Vector3i &value) {
	if (id == -1)
		return;
	glUniform3i(id, value.x, value.y, value.z);
}

void GLProgram::setParameter(int id, const Vector2 &value) {
	if (id == -1)
		return;
	glUniform2f(id, (GLfloat) value.x, (GLfloat) value.y);
}

void GLProgram::setParameter(int id, const Vector2i &value) {
	if (id == -1)
		return;
	glUniform2i(id, value.x, value.y);
}

void GLProgram::setParameter(int id, const Vector4 &value) {
	if (id == -1)
		return;
	glUniform4f(id, (GLfloat) value.x, (GLfloat) value.y, (GLfloat) value.z, (GLfloat) value.w);
}

void GLProgram::setParameter(int id, const GPUTexture *value) {
	if (id == -1)
		return;
	const std::set<int> &units = value->getTextureUnits();
	if (units.size() > 0) 
		glUniform1i(id, *(units.begin()));
	else
		Log(EWarn, "Unable to supply unbound texture \"%s\" to shader \"%s\"",
			value->getName().c_str(), getName().c_str());
}

void GLProgram::setParameter(int id, const Transform &trafo) {
	if (id == -1)
		return;
#ifdef SINGLE_PRECISION
	glUniformMatrix4fv(id, 1, true, reinterpret_cast<const GLfloat *>
		(trafo.getMatrix()->m));
#else
	GLfloat tmp[16];
	int i, j, idx=0;
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			tmp[idx++] = (GLfloat) trafo.getMatrix()->m[i][j];
	glUniformMatrix4fv(id, 1, true, tmp);
#endif
}

void GLProgram::setParameter(int id, const Spectrum &value) {
	if (id == -1)
		return;

	Float r, g, b;
	value.toLinearRGB(r, g, b);
	glUniform3f(id, r, g, b);
}

void GLProgram::bind() {
	glUseProgramObjectARB(m_program);
	m_bound = true;
}

void GLProgram::unbind() {
	glUseProgramObjectARB(0);
	m_bound = false;
}

void GLProgram::cleanup() {
	Log(EDebug, "Freeing GPU program \"%s\"", m_name.c_str());
	if (m_id[EVertexProgram] != 0) {
		if (m_program != 0)
			glDetachObjectARB(m_program, m_id[EVertexProgram]);
		glDeleteObjectARB(m_id[EVertexProgram]);
		m_id[EVertexProgram] = 0;
	}
	if (m_id[EFragmentProgram] != 0) {
		if (m_program != 0)
			glDetachObjectARB(m_program, m_id[EFragmentProgram]);
		glDeleteObjectARB(m_id[EFragmentProgram]);
		m_id[EFragmentProgram] = 0;
	}
	if (m_id[EGeometryProgram] != 0) {
		if (m_program != 0)
			glDetachObjectARB(m_program, m_id[EGeometryProgram]);
		glDeleteObjectARB(m_id[EGeometryProgram]);
		m_id[EGeometryProgram] = 0;
	}
	if (m_program != 0) {
		glDeleteObjectARB(m_program);
		m_program = 0;
	}
}

GLProgram::~GLProgram() {
	if (m_program != 0)
		cleanup();
}

MTS_IMPLEMENT_CLASS(GLProgram, false, GPUProgram)
MTS_NAMESPACE_END
