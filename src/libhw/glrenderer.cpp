/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#include <Carbon/Carbon.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gltexture.h>
#include <mitsuba/hw/glgeometry.h>
#include <mitsuba/hw/glprogram.h>
#include <mitsuba/hw/glsync.h>
#include <mitsuba/hw/font.h>

static mitsuba::PrimitiveThreadLocal<GLEWContextStruct> glewContext;

GLEWContextStruct *glewGetContext() {
	return &glewContext.get();
}

MTS_NAMESPACE_BEGIN

GLRenderer::GLRenderer(Session *session)
 : Renderer(session) {
}

GLRenderer::~GLRenderer() {
}

void GLRenderer::init(Device *device, Renderer *other) {
	Renderer::init(device, other);

	m_driverRenderer = (char *) glGetString(GL_RENDERER);
	m_driverVendor = (char *) glGetString(GL_VENDOR);
	m_driverVersion = (char *) glGetString(GL_VERSION);

	Log(m_logLevel, "OpenGL renderer : %s", m_driverRenderer.c_str());
	Log(m_logLevel, "OpenGL vendor   : %s", m_driverVendor.c_str());
	Log(m_logLevel, "OpenGL version  : %s", m_driverVersion.c_str());

	/* OpenGL extensions */
	GLenum err = glewInit();
	if (err != GLEW_OK) 
		Log(EError, "GLEW Error: %s\n", glewGetErrorString(err));

	if (glewIsSupported("GL_EXT_framebuffer_object")) {
		m_capabilities->setSupported(
			RendererCapabilities::ERenderToTexture, true);
		Log(m_logLevel, "Capabilities: Framebuffers objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Framebuffers objects are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_shading_language_100")) {
		m_capabilities->setSupported(
			RendererCapabilities::EShadingLanguage, true);
		Log(m_logLevel, "Capabilities: GLSL is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: GLSL is NOT supported!");
	}

	if (glewIsSupported("GL_ARB_texture_float")) {
		m_capabilities->setSupported(
			RendererCapabilities::EFloatingPointTextures, true);
		Log(m_logLevel, "Capabilities: Floating point textures are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Floating point textures are NOT supported!");
	}

	bool leopardWorkaround = false;
#if defined(__OSX__)
	/* Floating point color render buffers sort-of work for
	   Leopard/8600M or 9600M, but the extension is not reported */
	SInt32 MacVersion;
	if (Gestalt(gestaltSystemVersion, &MacVersion) == noErr) {
		if (MacVersion >= 0x1050 && MacVersion < 0x1060) {
			Log(EInfo, "Enabling Leopard floating point color buffer workaround");
			leopardWorkaround = true;
		}
	}
#endif

	if (glewIsSupported("GL_ARB_color_buffer_float") || leopardWorkaround) {
		m_capabilities->setSupported(
			RendererCapabilities::EFloatingPointBuffer, true);
		Log(m_logLevel, "Capabilities: Floating point color buffers are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Floating point color buffers are NOT supported!");
	}

	if (glewIsSupported("GL_EXT_framebuffer_blit")) {
		m_capabilities->setSupported(
			RendererCapabilities::EBufferBlit, true);
		Log(m_logLevel, "Capabilities: Fast buffer blitting is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Fast buffer blitting is NOT supported!");
	}

	if (glewIsSupported("GL_EXT_framebuffer_multisample") && 
		glewIsSupported("GL_EXT_framebuffer_blit") &&
		glewIsSupported("GL_ARB_texture_multisample")) {
		m_capabilities->setSupported(
			RendererCapabilities::EMultisampleRenderToTexture, true);
		Log(m_logLevel, "Capabilities: Multisample framebuffer objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Multisample framebuffer objects are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_vertex_buffer_object")) {
		m_capabilities->setSupported(
			RendererCapabilities::EVertexBufferObjects, true);
		Log(m_logLevel, "Capabilities: Vertex buffer objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Vertex buffer objects are NOT supported!");
	}

	if (glewIsSupported("GL_EXT_geometry_shader4")) {
		m_capabilities->setSupported(
			RendererCapabilities::EGeometryShaders, true);
		Log(m_logLevel, "Capabilities: Geometry shaders are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Geometry shaders are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_sync")) {
		m_capabilities->setSupported(
			RendererCapabilities::ESyncObjects, true);
		Log(m_logLevel, "Capabilities: Synchronization objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Synchronization objects are NOT supported!");
	}

	if (glewIsSupported("GL_NV_vertex_buffer_unified_memory")) {
		m_capabilities->setSupported(
			RendererCapabilities::EBindless, true);
		Log(m_logLevel, "Capabilities: Bindless rendering is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Bindless rendering is NOT supported!");
	}

	/* Hinting */
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	/* Disable color value clamping */
	if (m_capabilities->isSupported(
			RendererCapabilities::EFloatingPointBuffer) && !leopardWorkaround) {
		glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
	}

	/* Clip to viewport */
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	setBlendMode(EBlendNone);

	glEnable(GL_POINT_SMOOTH);

	m_normalsEnabled = false;
	m_texcoordsEnabled = false;
	m_tangentsEnabled = false;
	m_colorsEnabled = false;
	m_stride = -1;
	m_queuedTriangles = 0;

	checkError();
}

void GLRenderer::shutdown() {
	Renderer::shutdown();
}

GPUTexture *GLRenderer::createGPUTexture(const std::string &name,
		Bitmap *bitmap) {
	return new GLTexture(name, bitmap);
}

GPUGeometry *GLRenderer::createGPUGeometry(const TriMesh *mesh) {
	return new GLGeometry(mesh);
}

GPUProgram *GLRenderer::createGPUProgram(const std::string &name) {
	return new GLProgram(name);
}
	
GPUSync *GLRenderer::createGPUSync() {
	return new GLSync();
}

void GLRenderer::clear() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	checkError();
}

void GLRenderer::checkError(bool onlyWarn) {
	int glError = glGetError();

	if (glError)
		Log(onlyWarn ? m_warnLogLevel : EError, "OpenGL Error : %s", gluErrorString(glError));
}

void GLRenderer::beginDrawingMeshes(bool transmitOnlyPositions) {
	m_transmitOnlyPositions = transmitOnlyPositions;
	glEnableClientState(GL_VERTEX_ARRAY);
	m_stride = -1;

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		glEnableClientState(GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV);
		glEnableClientState(GL_ELEMENT_ARRAY_UNIFIED_NV);
	}
}

void GLRenderer::drawTriMesh(const TriMesh *mesh) {
	std::map<const TriMesh *, GPUGeometry *>::iterator it = m_geometry.find(mesh);
	if (it != m_geometry.end()) {
		/* Draw using vertex buffer objects (bindless if supported) */
		GLGeometry *geometry = static_cast<GLGeometry *>((*it).second);
		if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
			glBufferAddressRangeNV(GL_VERTEX_ARRAY_ADDRESS_NV, 0, 
				geometry->m_vertexAddr, geometry->m_vertexSize);
			int stride = geometry->m_stride;
			if (stride != m_stride) {
				glVertexFormatNV(3, GL_FLOAT, stride);
				glNormalFormatNV(GL_FLOAT, stride);
				glClientActiveTexture(GL_TEXTURE0);
				glTexCoordFormatNV(2, GL_FLOAT, stride);
				glClientActiveTexture(GL_TEXTURE1);
				glTexCoordFormatNV(3, GL_FLOAT, stride);
				glColorFormatNV(3, GL_FLOAT, stride);
				m_stride = stride;
			}

			if (!m_transmitOnlyPositions) {
				int pos = 3 * sizeof(GLfloat);
				
				if (mesh->hasVertexNormals()) {
					if (!m_normalsEnabled) {
						glEnableClientState(GL_NORMAL_ARRAY);
						m_normalsEnabled = true;
					}
					glBufferAddressRangeNV(GL_NORMAL_ARRAY_ADDRESS_NV, 0, 
						geometry->m_vertexAddr + pos, 
						geometry->m_vertexSize - pos);

					pos += 3 * sizeof(GLfloat);
				} else if (!m_normalsEnabled) {
					glDisableClientState(GL_NORMAL_ARRAY);
					m_normalsEnabled = false;
				}

				if (mesh->hasVertexTexcoords()) {
					glClientActiveTexture(GL_TEXTURE0);
					if (!m_texcoordsEnabled) {
						glEnableClientState(GL_TEXTURE_COORD_ARRAY);
						m_texcoordsEnabled = true;
					}
					glBufferAddressRangeNV(GL_TEXTURE_COORD_ARRAY_ADDRESS_NV, 0,
						geometry->m_vertexAddr + pos,
						geometry->m_vertexSize - pos);

					pos += 2 * sizeof(GLfloat);
				} else if (m_texcoordsEnabled) {
					glClientActiveTexture(GL_TEXTURE0);
					glDisableClientState(GL_TEXTURE_COORD_ARRAY);
					m_texcoordsEnabled = false;
				}

				/* Pass 'dpdu' as second set of texture coordinates */
				if (mesh->hasVertexTangents()) {
					glClientActiveTexture(GL_TEXTURE1);
					if (!m_tangentsEnabled) {
						glEnableClientState(GL_TEXTURE_COORD_ARRAY);
						m_tangentsEnabled = true;
					}

					glBufferAddressRangeNV(GL_TEXTURE_COORD_ARRAY_ADDRESS_NV, 1,
						geometry->m_vertexAddr + pos,
						geometry->m_vertexSize - pos);
					pos += 3 * sizeof(GLfloat);
				} else if (m_tangentsEnabled) {
					glClientActiveTexture(GL_TEXTURE1);
					glDisableClientState(GL_TEXTURE_COORD_ARRAY);
					m_tangentsEnabled = false;
				}

				if (mesh->hasVertexColors()) {
					if (!m_colorsEnabled) {
						glEnableClientState(GL_COLOR_ARRAY);
						m_colorsEnabled = true;
					}

					glBufferAddressRangeNV(GL_COLOR_ARRAY_ADDRESS_NV, 0, 
						geometry->m_vertexAddr + pos,
						geometry->m_vertexSize - pos);
				} else if (m_colorsEnabled) {
					glDisableClientState(GL_COLOR_ARRAY);
					m_colorsEnabled = false;
				}
			}
			glBufferAddressRangeNV(GL_ELEMENT_ARRAY_ADDRESS_NV, 0, 
				geometry->m_indexAddr, geometry->m_indexSize);
		} else {
			glBindBuffer(GL_ARRAY_BUFFER, geometry->m_vertexID);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->m_indexID);
			int stride = geometry->m_stride;

			/* Set up the vertex/normal arrays */
			glVertexPointer(3, GL_FLOAT, stride, (GLfloat *) 0);

			if (!m_transmitOnlyPositions) {
				int pos = 3;
				if (mesh->hasVertexNormals()) {
					if (!m_normalsEnabled) {
						glEnableClientState(GL_NORMAL_ARRAY);
						m_normalsEnabled = true;
					}
					glNormalPointer(GL_FLOAT, stride, (GLfloat *) 0 + pos);
					pos += 3;
				} else if (m_normalsEnabled) {
					glDisableClientState(GL_NORMAL_ARRAY);
					m_normalsEnabled = false;
				}

				if (mesh->hasVertexTexcoords()) {
					glClientActiveTexture(GL_TEXTURE0);
					if (!m_texcoordsEnabled) {
						glEnableClientState(GL_TEXTURE_COORD_ARRAY);
						m_texcoordsEnabled = true;
					}
					glTexCoordPointer(2, GL_FLOAT, stride, (GLfloat *) 0 + pos);
					pos += 2;
				} else if (m_texcoordsEnabled) {
					glClientActiveTexture(GL_TEXTURE0);
					glDisableClientState(GL_TEXTURE_COORD_ARRAY);
					m_texcoordsEnabled = false;
				}

				/* Pass 'dpdu' as second set of texture coordinates */
				if (mesh->hasVertexTangents()) {
					glClientActiveTexture(GL_TEXTURE1);
					if (!m_tangentsEnabled) {
						glEnableClientState(GL_TEXTURE_COORD_ARRAY);
						m_tangentsEnabled = true;
					}
					glTexCoordPointer(3, GL_FLOAT, stride, (GLfloat *) 0 + pos);
					pos += 3;
				} else if (m_tangentsEnabled) {
					glClientActiveTexture(GL_TEXTURE1);
					glDisableClientState(GL_TEXTURE_COORD_ARRAY);
					m_tangentsEnabled = false;
				}

				if (mesh->hasVertexColors()) {
					if (!m_colorsEnabled) {
						glEnableClientState(GL_COLOR_ARRAY);
						m_colorsEnabled = true;
					}
					glColorPointer(3, GL_FLOAT, stride, (GLfloat *) 0 + pos);
				} else if (m_colorsEnabled) {
					glDisableClientState(GL_COLOR_ARRAY);
					m_colorsEnabled = false;
				}
			}
		}

		size_t size = mesh->getTriangleCount();
		if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
			/* Draw all triangles */
			glDrawElements(GL_TRIANGLES, (GLsizei) (size * 3), 
				GL_UNSIGNED_INT, (GLvoid *) 0);
			m_queuedTriangles += size;
		} else {
			/* Spoon-feed them (keeps the OS responsive) */
			size_t size = mesh->getTriangleCount(), cur = 0;
			while (cur < size) {
				size_t drawAmt = std::min(size - cur,
						MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
				if (drawAmt > 0)
					glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3), 
						GL_UNSIGNED_INT, (GLuint *) 0 + cur * 3);
				m_queuedTriangles += drawAmt; cur += drawAmt;
				if (cur < size) {
					finish();
				}
			}
		}

		if (!m_capabilities->isSupported(RendererCapabilities::EBindless)) {
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	} else {
		/* Draw the old-fashioned way without VBOs */
		const GLchar *positions = (const GLchar *) mesh->getVertexPositions();
		const GLchar *normals = (const GLchar *) mesh->getVertexNormals();
		const GLchar *texcoords = (const GLchar *) mesh->getVertexTexcoords();
		const GLchar *tangents = (const GLchar *) mesh->getVertexTangents();
		const GLchar *colors = (const GLchar *) mesh->getVertexColors();
		const GLchar *indices  = (const GLchar *) mesh->getTriangles();
		GLenum dataType = sizeof(Float) == 4 ? GL_FLOAT : GL_DOUBLE;

		glVertexPointer(3, dataType, 0, positions);

		if (!m_transmitOnlyPositions) {
			if (mesh->hasVertexNormals()) {
				if (!m_normalsEnabled) {
					glEnableClientState(GL_NORMAL_ARRAY);
					m_normalsEnabled = true;
				}
				glNormalPointer(dataType, 0, normals);
			} else if (m_normalsEnabled) {
				glDisableClientState(GL_NORMAL_ARRAY);
				m_normalsEnabled = false;
			}

			glClientActiveTexture(GL_TEXTURE0);
			if (mesh->hasVertexTexcoords()) {
				if (!m_texcoordsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_texcoordsEnabled = true;
				}
				glTexCoordPointer(2, dataType, 0, texcoords);
			} else if (m_texcoordsEnabled) {
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_texcoordsEnabled = false;
			}

			/* Pass 'dpdu' as second set of texture coordinates */
			glClientActiveTexture(GL_TEXTURE1);
			if (mesh->hasVertexTangents()) {
				if (!m_tangentsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_tangentsEnabled = true;
				}
				glTexCoordPointer(3, dataType, sizeof(Vector), tangents); 
			} else if (m_tangentsEnabled) {
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_tangentsEnabled = false;
			}

			if (mesh->hasVertexColors()) {
				if (!m_colorsEnabled) {
					glEnableClientState(GL_COLOR_ARRAY);
					m_colorsEnabled = true;
				}
				// This won't work for spectral rendering
				glColorPointer(3, dataType, 0, colors);
			} else if (m_colorsEnabled) {
				glDisableClientState(GL_COLOR_ARRAY);
				m_colorsEnabled = false;
			}
		}

		size_t size = mesh->getTriangleCount();
		if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
			/* Draw all triangles */
			glDrawElements(GL_TRIANGLES, (GLsizei) (mesh->getTriangleCount()*3), 
				GL_UNSIGNED_INT, indices);
			m_queuedTriangles += size;
		} else {
			/* Spoon-feed them (keeps the OS responsive) */
			size_t size = mesh->getTriangleCount(), cur = 0;
			while (cur < size) {
				size_t drawAmt = std::min(size - cur,
						MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
				if (drawAmt > 0)
					glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3), 
						GL_UNSIGNED_INT, indices + cur * 3);
				m_queuedTriangles += drawAmt; cur += drawAmt;
				if (cur < size) {
					finish();
				}
			}
		}
	}
}

void GLRenderer::endDrawingMeshes() {
	glDisableClientState(GL_VERTEX_ARRAY);
	if (m_normalsEnabled) {
		glDisableClientState(GL_NORMAL_ARRAY);
		m_normalsEnabled = false;
	}
	if (m_texcoordsEnabled) {
		glClientActiveTexture(GL_TEXTURE0);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		m_texcoordsEnabled = false;
	}
	if (m_tangentsEnabled) {
		glClientActiveTexture(GL_TEXTURE1);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		m_tangentsEnabled = false;
	}
	if (m_colorsEnabled) {
		glDisableClientState(GL_COLOR_ARRAY);
		m_colorsEnabled = false;
	}

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		glDisableClientState(GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV);
		glDisableClientState(GL_ELEMENT_ARRAY_UNIFIED_NV);
	} else {
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
}
	
void GLRenderer::drawAll() {
	GLRenderer::beginDrawingMeshes(true);
	std::map<const TriMesh *, GPUGeometry *>::iterator it;
	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		for (it = m_geometry.begin(); it != m_geometry.end(); ++it) {
			const TriMesh *mesh = static_cast<const TriMesh *>((*it).first);
			const GLGeometry *geometry = static_cast<const GLGeometry *>((*it).second);
			int stride = geometry->m_stride;
			if (stride != m_stride) {
				glVertexFormatNV(3, GL_FLOAT, stride);
				m_stride = stride;
			}

			glBufferAddressRangeNV(GL_VERTEX_ARRAY_ADDRESS_NV, 0, 
				geometry->m_vertexAddr, geometry->m_vertexSize);
			glBufferAddressRangeNV(GL_ELEMENT_ARRAY_ADDRESS_NV, 0, 
				geometry->m_indexAddr, geometry->m_indexSize);

			size_t size = mesh->getTriangleCount();

			if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
				/* Draw all triangles */
				glDrawElements(GL_TRIANGLES, (GLsizei) (size * 3), 
					GL_UNSIGNED_INT, (GLvoid *) 0);
				m_queuedTriangles += size;
			} else {
				/* Spoon-feed them (keeps the OS responsive) */
				size_t size = mesh->getTriangleCount(), cur = 0;
				while (cur < size) {
					size_t drawAmt = std::min(size - cur,
							MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
					if (drawAmt > 0)
						glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3), 
							GL_UNSIGNED_INT, (GLuint *) 0 + cur * 3);
					m_queuedTriangles += drawAmt; cur += drawAmt;
					if (cur < size)
						finish();
				}
			}
		}
	} else {
		for (it = m_geometry.begin(); it != m_geometry.end(); ++it) {
			const TriMesh *mesh = static_cast<const TriMesh *>((*it).first);
			const GLGeometry *geometry = static_cast<const GLGeometry *>((*it).second);

			glBindBuffer(GL_ARRAY_BUFFER, geometry->m_vertexID);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->m_indexID);

			/* Set up the vertex/normal arrays */
			glVertexPointer(3, GL_FLOAT, geometry->m_stride, (GLfloat *) 0);

			size_t size = mesh->getTriangleCount();

			if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
				/* Draw all triangles */
				glDrawElements(GL_TRIANGLES, (GLsizei) (size * 3), 
					GL_UNSIGNED_INT, (GLvoid *) 0);
				m_queuedTriangles += size;
			} else {
				/* Spoon-feed them (keeps the OS responsive) */
				size_t size = mesh->getTriangleCount(), cur = 0;
				while (cur < size) {
					size_t drawAmt = std::min(size - cur,
							MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
					if (drawAmt > 0)
						glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3), 
							GL_UNSIGNED_INT, (GLuint *) 0 + cur * 3);
					m_queuedTriangles += drawAmt; cur += drawAmt;
					if (cur < size)
						finish();
				}
			}
		}
	}
	GLRenderer::endDrawingMeshes();
}

void GLRenderer::blitTexture(const GPUTexture *tex, bool flipVertically,
		bool centerHoriz, bool centerVert, const Vector2i &offset) {
	tex->bind();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	if (tex->getType() == GPUTexture::ETexture2D) {
		GLint viewport[4];	
		glGetIntegerv(GL_VIEWPORT, viewport);
		Vector2i scrSize = Vector2i(viewport[2], viewport[3]);
		Vector2i texSize = Vector2i(tex->getSize().x, tex->getSize().y);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glBegin(GL_QUADS);

		Vector2i upperLeft(0), lowerRight(0);
		if (centerHoriz)
			upperLeft.x = (scrSize.x - texSize.x)/2;
		if (centerVert)
			upperLeft.y = (scrSize.y - texSize.y)/2;
		upperLeft += offset;
		lowerRight = upperLeft + texSize;

		if (flipVertically)
			std::swap(upperLeft.y, lowerRight.y);

		const float zDepth = -1.0f; // just before the far plane
		glTexCoord2f(0.0f, 0.0f);
		glVertex3f(upperLeft.x, upperLeft.y, zDepth);
		glTexCoord2f(1.0f, 0.0f);
		glVertex3f(lowerRight.x, upperLeft.y, zDepth);
		glTexCoord2f(1.0f, 1.0f);
		glVertex3f(lowerRight.x, lowerRight.y, zDepth);
		glTexCoord2f(0.0f, 1.0f);
		glVertex3f(upperLeft.x, lowerRight.y, zDepth);
		glEnd();
	} else if (tex->getType() == GPUTexture::ETextureCubeMap) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		/* From OpenTK */
		glBegin(GL_QUADS);
		// 0 -x
		glTexCoord3f(-1.0f, +1.0f, -1.0f);
		glVertex2f(-1.0f, +0.333f);
		glTexCoord3f(-1.0f, +1.0f, +1.0f);
		glVertex2f(-0.5f, +0.333f);
		glTexCoord3f(-1.0f, -1.0f, +1.0f);
		glVertex2f(-0.5f, -0.333f);
		glTexCoord3f(-1.0f, -1.0f, -1.0f);
		glVertex2f(-1.0f, -0.333f);

		// 1 +z
		glTexCoord3f(-1.0f, +1.0f, +1.0f);
		glVertex2f(-0.5f, +0.333f);
		glTexCoord3f(+1.0f, +1.0f, +1.0f);
		glVertex2f(+0.0f, +0.333f);
		glTexCoord3f(+1.0f, -1.0f, +1.0f);
		glVertex2f(+0.0f, -0.333f);
		glTexCoord3f(-1.0f, -1.0f, +1.0f);
		glVertex2f(-0.5f, -0.333f);

		// 2 +x
		glTexCoord3f(+1.0f, +1.0f, +1.0f);
		glVertex2f(+0.0f, +0.333f);
		glTexCoord3f(+1.0f, +1.0f, -1.0f);
		glVertex2f(+0.5f, +0.333f);
		glTexCoord3f(+1.0f, -1.0f, -1.0f);
		glVertex2f(+0.5f, -0.333f);
		glTexCoord3f(+1.0f, -1.0f, +1.0f);
		glVertex2f(+0.0f, -0.333f);

		// 3 -z
		glTexCoord3f(+1.0f, +1.0f, -1.0f);
		glVertex2f(+0.5f, +0.333f);
		glTexCoord3f(-1.0f, +1.0f, -1.0f);
		glVertex2f(+1.0f, +0.333f);
		glTexCoord3f(-1.0f, -1.0f, -1.0f);
		glVertex2f(+1.0f, -0.333f);
		glTexCoord3f(+1.0f, -1.0f, -1.0f);
		glVertex2f(+0.5f, -0.333f);

		// 4 +y
		glTexCoord3f(-1.0f, +1.0f, -1.0f);
		glVertex2f(-0.5f, +1.0f);
		glTexCoord3f(+1.0f, +1.0, -1.0f);
		glVertex2f(+0.0f, +1.0);
		glTexCoord3f(+1.0f, +1.0f, +1.0f);
		glVertex2f(+0.0f, +0.333f);
		glTexCoord3f(-1.0f, +1.0f, +1.0f);
		glVertex2f(-0.5f, +0.333f);

		// 5 -y
		glTexCoord3f(-1.0f, -1.0f, +1.0f);
		glVertex2f(-0.5f, -0.333f);
		glTexCoord3f(+1.0f, -1.0, +1.0f);
		glVertex2f(+0.0f, -0.333f);
		glTexCoord3f(+1.0f, -1.0f, -1.0f);
		glVertex2f(+0.0f, -1.0f);
		glTexCoord3f(-1.0f, -1.0f, -1.0f);
		glVertex2f(-0.5f, -1.0f);
		glEnd();
	}
	tex->unbind();
}

void GLRenderer::blitQuad(bool flipVertically) {
	GLint viewport[4];	
	glGetIntegerv(GL_VIEWPORT, viewport);
	Vector2 scrSize(viewport[2], viewport[3]);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	const Float zDepth = -1.0f;
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, flipVertically ? 1.0f : 0.0f);
	glVertex3f(0.0f, 0.0f, zDepth);
	glTexCoord2f(1.0f, flipVertically ? 1.0f : 0.0f);
	glVertex3f(scrSize.x, 0.0f, zDepth);
	glTexCoord2f(1.0f, flipVertically ? 0.0f : 1.0f);
	glVertex3f(scrSize.x, scrSize.y, zDepth);
	glTexCoord2f(0.0f, flipVertically ? 0.0f : 1.0f);
	glVertex3f(0.0f, scrSize.y, zDepth);
	glEnd();
}

void GLRenderer::drawText(const Point2i &_pos, 
		const Font *font, const std::string &text) {
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	Vector2i scrSize = Vector2i(viewport[2], viewport[3]);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	font->getTexture()->bind();
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	Point2i pos(_pos);
	int initial = pos.x;

	glBegin(GL_QUADS);
	for (size_t i=0; i<text.length(); i++) {
		char character = text[i];
		if (character == '\r')
			continue;
		if (character == '\n') {
			pos.x = initial;
			pos.y += (int) (font->getMaxVerticalBearing()*4.0/3.0);
			continue;
		}

		const Font::Glyph &glyph = font->getGlyph(character);

		Point2i start = pos + Vector2i(
			glyph.horizontalBearing,
			font->getMaxVerticalBearing() - glyph.verticalBearing
		);
		Point2i end = start + glyph.size;
		Point2 txStart = glyph.tx;
		Point2 txEnd = txStart + glyph.ts;

		glTexCoord2f(txStart.x, txStart.y);
		glVertex2f(    start.x,   start.y);
		glTexCoord2f(txEnd.x,   txStart.y);
		glVertex2f(    end.x,     start.y);
		glTexCoord2f(txEnd.x,     txEnd.y);
		glVertex2f(    end.x,       end.y);
		glTexCoord2f(txStart.x,   txEnd.y);
		glVertex2f(    start.x,     end.y);

		pos.x += glyph.horizontalAdvance;

		if (i+1 < text.length())
			pos.x += font->getKerning(character, text[i+1]);
	}
	glEnd();

	font->getTexture()->unbind();
	glDisable(GL_BLEND);
}

void GLRenderer::setPointSize(Float size) {
	glPointSize(size);
}

void GLRenderer::drawPoint(const Point &p) {
	glBegin(GL_POINTS);
	glVertex3f(p.x, p.y, p.z);
	glEnd();
}

void GLRenderer::drawLine(const Point &a, const Point &b) {
	glBegin(GL_LINES);
	glVertex3f(a.x, a.y, a.z);
	glVertex3f(b.x, b.y, b.z);
	glEnd();
}

void GLRenderer::drawEllipse(const Point &center, 
		const Vector &axis1, const Vector &axis2) {
	const int nSteps = 100;
	const float stepSize = 2*M_PI/nSteps;
	glBegin(GL_LINE_LOOP);
	for (int i=0; i<100; ++i) {
		Point p = center + axis1 * std::cos(i*stepSize) 
			+ axis2 * std::sin(i*stepSize);
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
}

void GLRenderer::drawAABB(const AABB &aabb) {
	#define V(a,b,c) glVertex3f(aabb.a.x, aabb.b.y, aabb.c.z)
	glBegin(GL_LINE_LOOP); V(max,min,max); V(max,min,min); V(max,max,min); V(max,max,max); glEnd();
	glBegin(GL_LINE_LOOP); V(max,max,max); V(max,max,min); V(min,max,min); V(min,max,max); glEnd();
	glBegin(GL_LINE_LOOP); V(max,max,max); V(min,max,max); V(min,min,max); V(max,min,max); glEnd();
	glBegin(GL_LINE_LOOP); V(min,min,max); V(min,max,max); V(min,max,min); V(min,min,min); glEnd();
	glBegin(GL_LINE_LOOP); V(min,min,max); V(min,min,min); V(max,min,min); V(max,min,max); glEnd();
	glBegin(GL_LINE_LOOP); V(min,min,min); V(min,max,min); V(max,max,min); V(max,min,min); glEnd();
	#undef V
}

void GLRenderer::setCamera(const ProjectiveCamera *camera) {
	GLfloat temp1[16], temp2[16];
	const Matrix4x4 &view = camera->getViewTransform().getMatrix();
	const Matrix4x4 &proj = camera->getGLProjectionTransform().getMatrix();

	int pos=0;
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			temp1[pos]=(GLfloat) proj.m[i][j];
			temp2[pos++]=(GLfloat) view.m[i][j];
		}
	}

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(temp1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0f, 1.0f, -1.0f);
	glMultMatrixf(temp2);
}

void GLRenderer::setCamera(const ProjectiveCamera *camera, const Point2 &jitter) {
	GLfloat temp1[16], temp2[16];
	const Matrix4x4 &view = camera->getViewTransform().getMatrix();
	const Matrix4x4 &proj = camera->getGLProjectionTransform(jitter).getMatrix();

	int pos=0;
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			temp1[pos]=(GLfloat) proj.m[i][j];
			temp2[pos++]=(GLfloat) view.m[i][j];
		}
	}

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(temp1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0f, 1.0f, -1.0f);
	glMultMatrixf(temp2);
}

void GLRenderer::setCamera(const Matrix4x4 &proj, const Matrix4x4 &view) {
	GLfloat temp1[16], temp2[16];
	int pos=0;
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			temp1[pos]=(GLfloat) proj.m[i][j];
			temp2[pos++]=(GLfloat) view.m[i][j];
		}
	}

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(temp1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0f, 1.0f, -1.0f);
	glMultMatrixf(temp2);
}

void GLRenderer::setDepthOffset(Float value) {
	if (value == 0)
		glDisable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(4.0f, (GLfloat) value);
}

void GLRenderer::setDepthMask(bool value) {
	glDepthMask(value ? GL_TRUE : GL_FALSE);
}

void GLRenderer::setDepthTest(bool value) {
	if (value)
		glEnable(GL_DEPTH_TEST);
	else
		glDisable(GL_DEPTH_TEST);
}

void GLRenderer::setColorMask(bool value) {
	GLboolean flag = value ? GL_TRUE : GL_FALSE;
	glColorMask(flag, flag, flag, flag);
}

void GLRenderer::flush() {
	glFlush();
}

void GLRenderer::finish() {
	glFinish();
	m_queuedTriangles = 0;
}

void GLRenderer::setColor(const Spectrum &spec) {
	Float r, g, b;
	spec.toLinearRGB(r, g, b);
	glColor4f(r, g, b, 1);
}

void GLRenderer::setBlendMode(EBlendMode mode) {
	switch (mode) {
		case EBlendNone:
			glDisable(GL_BLEND);
			break;
		case EBlendAdditive:
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
			break;
		case EBlendAlpha:
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			break;
		default:
			Log(EError, "Invalid blend mode!");
	}
}

void GLRenderer::setCullMode(ECullMode mode) {
	switch (mode) {
		case ECullNone:
			glDisable(GL_CULL_FACE);
			break;
		case ECullFront:
			glEnable(GL_CULL_FACE);
			glCullFace(GL_FRONT);
			break;
		case ECullBack:
			glEnable(GL_CULL_FACE);
			glCullFace(GL_BACK);
			break;
		default:
			Log(EError, "Invalid culling mode!");
	}
}

MTS_IMPLEMENT_CLASS(GLRenderer, true, Renderer)
MTS_NAMESPACE_END
