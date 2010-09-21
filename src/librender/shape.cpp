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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>

MTS_NAMESPACE_BEGIN

Shape::Shape(const Properties &props) 
 : ConfigurableObject(props), m_luminaire(NULL) {
	m_objectToWorld = props.getTransform("toWorld", Transform());
	m_worldToObject = m_objectToWorld.inverse();
	m_surfaceArea = 0.0f;
}

Shape::Shape(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_worldToObject = Transform(stream);
	m_aabb = AABB(stream);
	m_bsphere = BSphere(stream);
	m_surfaceArea = stream->readFloat();
	m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));
	m_subsurface = static_cast<Subsurface *>(manager->getInstance(stream));
	m_luminaire = static_cast<Luminaire *>(manager->getInstance(stream));
	m_invSurfaceArea = 1.0f / m_surfaceArea;
	m_objectToWorld = m_worldToObject.inverse();
}

Shape::~Shape() {
}


void Shape::configure() {
	/* Ensure that there is at least some default BSDF */
	if (m_bsdf == NULL) {
		m_bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(BSDF::m_theClass, Properties("lambertian")));
	}
}
	
bool Shape::isCompound() const {
	return false;
}

bool Shape::isClippable() const {
	return false;
}

AABB Shape::getClippedAABB(const AABB &aabb) const {
	AABB result(m_aabb);
	result.clip(aabb);
	return result;
}

Shape *Shape::getElement(int i) {
	return NULL;
}

Float Shape::pdfArea(const ShapeSamplingRecord &sRec) const {
	return m_invSurfaceArea;
}

Float Shape::sampleSolidAngle(ShapeSamplingRecord &sRec, const Point &from, const Point2 &sample) const {
	/* Turns the area sampling routine into one that samples wrt. solid angles */
	Float pdfArea = sampleArea(sRec, sample);
	Vector lumToPoint = from - sRec.p;
	Float distSquared = lumToPoint.lengthSquared(), dp = dot(lumToPoint, sRec.n);
	if (dp > 0) {
		return pdfArea * distSquared * std::sqrt(distSquared) / dp;
	} else {
		return 0.0f;
	}
}

Float Shape::pdfSolidAngle(const ShapeSamplingRecord &sRec, const Point &from) const {
	/* Turns the area sampling routine into one that samples wrt. solid angles */
	Vector lumToPoint = from - sRec.p;
	Float distSquared = lumToPoint.lengthSquared();
	Float invDP = std::max((Float) 0, std::sqrt(distSquared) / dot(lumToPoint, sRec.n));
	return pdfArea(sRec) * distSquared * invDP;
}

void Shape::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();
	if (cClass->derivesFrom(BSDF::m_theClass)) {
		m_bsdf = static_cast<BSDF *>(child);
	} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
		Assert(m_luminaire == NULL);
		m_luminaire = static_cast<Luminaire *>(child);
	} else if (cClass->derivesFrom(Subsurface::m_theClass)) {
		Assert(m_subsurface == NULL);
		m_subsurface = static_cast<Subsurface *>(child);
	} else {
		Log(EError, "Shape: Invalid child node!");
	}
}

/// Ray intersection test
bool Shape::rayIntersect(const Ray &ray, Intersection &its) const {
	Log(EError, "Not implemented!");
	return false;
}

bool Shape::rayIntersect(const Ray &ray, Float start, Float end, Float &t) const {
	Log(EError, "Not implemented!");
	return false;
}

Float Shape::sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
	Log(EError, "Not implemented!");
	return 0;
}

#if defined(MTS_SSE)
__m128 Shape::rayIntersectPacket(const RayPacket4 &packet, const
       __m128 mint, __m128 maxt, __m128 inactive, Intersection4 &its) const {
	Log(EError, "Not implemented!");
	return SSEConstants::zero.ps;
}
#endif

void Shape::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);

	m_worldToObject.serialize(stream);
	m_aabb.serialize(stream);
	m_bsphere.serialize(stream);
	stream->writeFloat(m_surfaceArea);
	manager->serialize(stream, m_bsdf.get());
	manager->serialize(stream, m_subsurface.get());
	manager->serialize(stream, m_luminaire.get());
}

std::string ShapeSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "ShapeSamplingRecord[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  n = " << n.toString() << std::endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Shape, true, ConfigurableObject)
MTS_NAMESPACE_END
