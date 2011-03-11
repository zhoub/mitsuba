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

#include <mitsuba/render/camera.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter cameraRays("General", "Camera ray generations");

/**
 * Perspective camera model. Depth of field can optionally 
 * be activated by specifying both focal distance and lens 
 * radius (uses a thin lens approximation).
 * Based on the PBRT camera implementation.
 */
class PerspectiveCamera : public PinholeCamera {
public:
	PerspectiveCamera(const Properties &props) 
		: PinholeCamera(props) {
		/* Distance to the focal plane */
		m_focalDistance = props.getFloat("focalDistance", m_farClip);
		/* World-space lens radius */
		m_lensRadius = props.getFloat("lensRadius", 0.0f);
		if (m_cameraToWorld.hasScale()) 
			Log(EError, "Mitsuba's perspective camera does not allow scale "
				"factors in the camera-to-world transformation! Please remove those factors, "
				"since they will cause inconsistencies in various parts of the renderer.");
	}

	PerspectiveCamera(Stream *stream, InstanceManager *manager) 
	 : PinholeCamera(stream, manager) {
		m_worldToScreen = Transform(stream);
		m_rasterToCamera = Transform(stream);
		m_cameraToRaster = Transform(stream);
		m_rasterToScreen = Transform(stream);
		m_screenToRaster = Transform(stream);
		m_lensRadius = stream->readFloat();
		m_focalDistance = stream->readFloat();
		configure();
	}

	void configure() {
		PinholeCamera::configure();

		bool mapYToNDC01 = (m_aspect >= 1.0f);
		if (!m_mapSmallerSide)
			mapYToNDC01 = !mapYToNDC01;
		
		/* Maps from the image plane space to raster space. The
		   smaller 2D coordinate axis on the film will have
		   normalized device coordinates in [0, 1] (unless 
		   the flag mapSmallerSide is specified) . Also inverts
		   the Y coordinate
		*/
		m_cameraToScreen = Transform::perspective(m_fov, m_nearClip, m_farClip);

		if (mapYToNDC01) {
			m_screenToRaster = 
				Transform::scale(Vector((Float) m_film->getSize().x, (Float) m_film->getSize().y, 1.0f))
				* Transform::scale(Vector(1/(2*m_aspect), -0.5f, 1.0f))
				* Transform::translate(Vector(m_aspect, -1.0f, 0));
		} else {
			m_screenToRaster = 
				Transform::scale(Vector((Float) m_film->getSize().x, (Float) m_film->getSize().y, 1.0f))
				* Transform::scale(Vector(0.5f, -0.5f * m_aspect, 1.0f))
				* Transform::translate(Vector(1.0f, - 1 / m_aspect, 0));
		}

		m_cameraToScreenGL = Transform::glPerspective(m_yfov, m_nearClip, m_farClip)
			* Transform::scale(Vector(1/m_aspect, 1.0f, 1.0f));

		m_rasterToScreen = m_screenToRaster.inverse();

		m_worldToScreen = m_cameraToScreen * m_worldToCamera;
		m_rasterToCamera = m_cameraToScreen.inverse() * m_rasterToScreen;
		m_cameraToRaster = m_rasterToCamera.inverse();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		PinholeCamera::serialize(stream, manager);

		m_worldToScreen.serialize(stream);
		m_rasterToCamera.serialize(stream);
		m_cameraToRaster.serialize(stream);
		m_rasterToScreen.serialize(stream);
		m_screenToRaster.serialize(stream);
		stream->writeFloat(m_lensRadius);
		stream->writeFloat(m_focalDistance);
	}

	bool needsLensSample() const {
		return m_lensRadius > 0.0f;
	}

	void generateRay(const Point2 &dirSample, const Point2 &lensSample,
			Float timeSample, Ray &ray) const {
		++cameraRays;

		/* Calculate intersection on the image plane */
		Point rasterCoords(dirSample.x, dirSample.y, 0);
		Point imageCoords;
		m_rasterToCamera(rasterCoords, imageCoords);

		/* Construct ray in camera space */
		Ray localRay(Point(0, 0, 0), Vector(imageCoords),
			m_shutterOpen + m_shutterOpenTime * timeSample);

		if (m_lensRadius > 0.0f) {
			/* Sample a point on the aperture */
			Point2 lensPos = squareToDiskConcentric(lensSample)
				* m_lensRadius;

			/* Calculate the intersection with the focal plane */
			Point itsFocal = 
				localRay(m_focalDistance / localRay.d.z);

			/* Perturb the ray accordingly */
			localRay.o.x += lensPos.x;
			localRay.o.y += lensPos.y;
			localRay.d = itsFocal - localRay.o;
		}

		localRay.d = normalize(localRay.d);
		Float invZ = 1.0f / localRay.d.z;
		localRay.mint = m_nearClip * invZ;
		localRay.maxt = m_farClip * invZ;

		/* Transform into world space */
		m_cameraToWorld(localRay, ray);
	}

	bool positionToSample(const Point &p, Point2 &sample) const {
		Point localP, imageP;
		m_worldToCamera(p, localP);
		if (localP.z <= 0.0f)
			return false;
		m_cameraToRaster(localP, imageP);
		sample.x = imageP.x;
		sample.y = imageP.y;
		const Point2i max(m_film->getCropOffset() + m_film->getCropSize());
		return imageP.z > 0 && imageP.z < 1
			&& sample.x >= m_film->getCropOffset().x 
			&& sample.y >= m_film->getCropOffset().y
			&& sample.x < max.x && sample.y < max.y;
	}

	Transform getGLProjectionTransform(const Point2 &jitter) const {
		Float top = std::tan(m_yfov * M_PI/360) * m_nearClip, bottom = -top;
		Float left = m_aspect * bottom, right = m_aspect * top;
		Vector2 scale(
			(right-left)/m_film->getSize().x,
			(top-bottom)/m_film->getSize().y);
		return Transform::glFrustum(
			left - jitter.x*scale.y, right - jitter.x*scale.x, 
			bottom - jitter.y*scale.y, top - jitter.y*scale.y, 
			m_nearClip, m_farClip);
	}

	Float areaDensity(const Point2 &p) const {
		return 0.0f;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PerspectiveCamera[" << std::endl
			<< "  sampler = " << indent(m_sampler->toString()) << "," << std::endl
			<< "  film = " << indent(m_film->toString()) << "," << std::endl
			<< "  aspect = " << m_aspect << "," << std::endl
			<< "  xfov = " << m_xfov << "," << std::endl
			<< "  yfov = " << m_yfov << "," << std::endl
			<< "  nearClip = " << m_nearClip << "," << std::endl
			<< "  farClip = " << m_farClip << "," << std::endl
			<< "  shutterOpen = " << m_shutterOpen << "," << std::endl
			<< "  shutterClose = " << m_shutterClose << "," << std::endl
			<< "  lensRadius = " << m_lensRadius << "," << std::endl
			<< "  focalDistance = " << m_focalDistance << "," << std::endl
			<< "  cameraToWorld = " << indent(m_cameraToWorld.toString()) << "," << std::endl
			<< "  cameraToScreen = " << indent(m_cameraToScreen.toString()) << "," << std::endl
			<< "  rasterToCamera = " << indent(m_rasterToCamera.toString()) << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_worldToScreen;
	Transform m_rasterToCamera, m_cameraToRaster;
	Transform m_rasterToScreen, m_screenToRaster;
	Float m_lensRadius, m_focalDistance;
};


MTS_IMPLEMENT_CLASS_S(PerspectiveCamera, false, PinholeCamera)
MTS_EXPORT_PLUGIN(PerspectiveCamera, "Perspective camera");
MTS_NAMESPACE_END
