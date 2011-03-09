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
#include <mitsuba/render/medium.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

Camera::Camera(const Properties &props)
 : ConfigurableObject(props), m_properties(props) {
	m_cameraToWorld = props.getTransform("toWorld", Transform());
	m_shutterOpen = props.getFloat("shutterOpen", 0.0f);
	m_shutterClose = props.getFloat("shutterClose", 0.0f);
	if (m_shutterOpen > m_shutterClose)
		Log(EError, "Shutter opening time must be less than "
			"or equal to the shutter closing time!");
	m_worldToCamera = m_cameraToWorld.inverse();
	m_position = m_cameraToWorld(Point(0,0,0));
	m_shutterOpenTime = m_shutterClose - m_shutterOpen;
}

Camera::Camera(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
	m_film = static_cast<Film *>(manager->getInstance(stream));
	m_sampler = static_cast<Sampler *>(manager->getInstance(stream));
	m_medium = static_cast<Medium *>(manager->getInstance(stream));
	m_worldToCamera = Transform(stream);
	m_cameraToWorld = Transform(stream);
	m_shutterOpen = stream->readFloat();
	m_shutterClose = stream->readFloat();
	m_position = m_cameraToWorld(Point(0,0,0));
	m_shutterOpenTime = m_shutterClose - m_shutterOpen;
}

Camera::~Camera() {
}

void Camera::setParent(ConfigurableObject *parent) {
	// Do not store the parent - this is useful in case
	// the camera subtree needs to be serialized by itself.
}

Point Camera::getPosition(const Point2 &sample) const {
	return m_position; // default impl.
}

void Camera::setFilm(Film *film) {
	m_film = film;
	configure();
}

void Camera::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	manager->serialize(stream, m_film.get());
	manager->serialize(stream, m_sampler.get());
	manager->serialize(stream, m_medium.get());
	m_worldToCamera.serialize(stream);
	m_cameraToWorld.serialize(stream);
	stream->writeFloat(m_shutterOpen);
	stream->writeFloat(m_shutterClose);
}

void Camera::generateRayDifferential(const Point2 &sample,
	const Point2 &lensSample, Float timeSample, RayDifferential &ray) const {

	generateRay(sample, lensSample, timeSample, ray);
	Point2 temp = sample; temp.x += 1; 
	generateRay(temp, lensSample, timeSample, ray.rx);
	temp = sample; temp.y += 1;
	generateRay(temp, lensSample, timeSample, ray.ry);
	ray.hasDifferentials = true;
}

void Camera::addChild(const std::string &name, ConfigurableObject *child) {
	if (child->getClass()->derivesFrom(MTS_CLASS(Sampler))) {
		Assert(m_sampler == NULL);
		m_sampler = static_cast<Sampler *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Film))) {
		Assert(m_film == NULL);
		m_film = static_cast<Film *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Medium))) {
		Assert(m_medium == NULL);
		m_medium = static_cast<Medium *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Medium))) {
		Assert(m_medium == NULL);
		m_medium = static_cast<Medium *>(child);
	} else {
		Log(EError, "Camera: Invalid child node!");
	}
}

ProjectiveCamera::ProjectiveCamera(Stream *stream, InstanceManager *manager)
 : Camera(stream, manager) {
	m_cameraToScreen = Transform(stream);
	m_cameraToScreenGL = Transform(stream);
	m_aspect = stream->readFloat();
}
	
void ProjectiveCamera::configure() {
	if (m_film == NULL) {
		/* Instantiate an EXR film by default */
		m_film = static_cast<Film*> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Film), Properties("exrfilm")));
		m_film->configure();
	}

	if (m_sampler == NULL) {
		/* No sampler has been selected - load an independent filter with 4 samples/pixel by default */
		Properties props("independent");
		props.setInteger("sampleCount", 4);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();
	}
	m_aspect = (Float) m_film->getSize().x / (Float) m_film->getSize().y;
}

void ProjectiveCamera::serialize(Stream *stream, InstanceManager *manager) const {
	Camera::serialize(stream, manager);

	m_cameraToScreen.serialize(stream);
	m_cameraToScreenGL.serialize(stream);
	stream->writeFloat(m_aspect);
}

PinholeCamera::PinholeCamera(const Properties &props)
 : ProjectiveCamera(props) {
	/* Field of view of the camera (in degrees) */
	m_fov = props.getFloat("fov", 90);
	/* Specifies which side of the image plane should cover the
	   field of view specified in the <tt>fov</tt> parameter
	*/
	m_mapSmallerSide = props.getBoolean("mapSmallerSide", true);
}
	
PinholeCamera::PinholeCamera(Stream *stream, InstanceManager *manager)
 : ProjectiveCamera(stream, manager) {
	m_fov = stream->readFloat();
	m_mapSmallerSide = stream->readBool();
}

void PinholeCamera::serialize(Stream *stream, InstanceManager *manager) const {
	ProjectiveCamera::serialize(stream, manager);

	stream->writeFloat(m_fov);
	stream->writeBool(m_mapSmallerSide);
}

void PinholeCamera::configure() {
	ProjectiveCamera::configure();

	bool mapYToNDC01 = (m_aspect >= 1.0f);
	if (!m_mapSmallerSide)
		mapYToNDC01 = !mapYToNDC01;

	if (mapYToNDC01) {
		m_yfov = m_fov;
		m_xfov = radToDeg(2 * std::atan(std::tan(degToRad(m_fov)/2) * m_aspect));
	} else {
		m_xfov = m_fov;
		m_yfov = radToDeg(2 * std::atan(std::tan(degToRad(m_fov)/2) / m_aspect));
	}

	m_imagePlaneSize.x = 2 * std::tan(degToRad(m_xfov)/2);
	m_imagePlaneSize.y = 2 * std::tan(degToRad(m_yfov)/2);
	m_imagePlanePixelSize.x = m_imagePlaneSize.x / getFilm()->getSize().x;
	m_imagePlanePixelSize.y = m_imagePlaneSize.y / getFilm()->getSize().y;
	m_imagePlaneInvArea = 1 / (m_imagePlaneSize.x * m_imagePlaneSize.y);
}

Float PinholeCamera::importance(const Point2 &p) const {
	Float x = (p.x * m_imagePlanePixelSize.x) - .5f * m_imagePlaneSize.x;
	Float y = (p.y * m_imagePlanePixelSize.y) - .5f * m_imagePlaneSize.y;

	return std::pow(1 + x*x+y*y, (Float) (3.0/2.0)) * m_imagePlaneInvArea;
}

Float PinholeCamera::importance(const Vector &v) const {
	Vector localV;
	m_worldToCamera(v, localV);
	if (localV.z <= 0.0f) 
		return 0.0f;
	Float invZ = 1.0f / localV.z;
	localV.x *= invZ; localV.y *= invZ;
	if (2*std::abs(localV.x)>m_imagePlaneSize.x 
	 || 2*std::abs(localV.y)>m_imagePlaneSize.y) 
		return 0.0f;
	return std::pow(1 + localV.x*localV.x+localV.y*localV.y, 
		(Float) (3.0/2.0)) * m_imagePlaneInvArea;
}

Vector2 PinholeCamera::getImagePlaneSize(Float dist) const {
	return m_imagePlaneSize * dist;
}

MTS_IMPLEMENT_CLASS(Camera, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(ProjectiveCamera, true, Camera)
MTS_IMPLEMENT_CLASS(PinholeCamera, true, ProjectiveCamera)
MTS_NAMESPACE_END
