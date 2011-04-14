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

#if !defined(__GLXRENDERER_H)
#define __GLXRENDERER_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/glxdevice.h>

MTS_NAMESPACE_BEGIN

/** \brief GLX (XFree86) renderer
 */
class MTS_EXPORT_HW GLXRenderer : public GLRenderer {
public:
	/// Create a new renderer
	GLXRenderer(X11Session *session);

	/// Return the rendering context
	inline GLXContext getGLXContext() { return m_context; }
	
	/// Initialize the renderer
	void init(Device *device, Renderer *other = NULL);

	/// Shut the renderer down
	void shutdown();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLXRenderer();
private:
    GLXContext m_context;
};

MTS_NAMESPACE_END

#endif /* __GLXRENDERER_H */
