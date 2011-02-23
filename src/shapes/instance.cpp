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

#include "shapegroup.h"

MTS_NAMESPACE_BEGIN

class Instance : public Shape {
public:
	Instance(const Properties &props) : Shape(props) {
		m_objectToWorld = props.getTransform("toWorld", Transform());
		m_worldToObject = m_objectToWorld.inverse();
	}

	Instance(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		m_shapeGroup = static_cast<ShapeGroup *>(manager->getInstance(stream));
		m_objectToWorld = Transform(stream);
		m_worldToObject = m_objectToWorld.inverse();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
		manager->serialize(stream, m_shapeGroup.get());
	}

	void configure() {
		if (!m_shapeGroup)
			Log(EError, "A reference to a 'shapegroup' must be specified!");
	}

	BoundingBox3 getBoundingBox3() const {
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		const BoundingBox3 &aabb = kdtree->getBoundingBox3();
		if (!aabb.isValid()) // the geometry group is empty
			return aabb;
		BoundingBox3 result;
		for (int i=0; i<8; ++i)
			result.expandBy(m_objectToWorld(aabb.getCorner(i)));
		return result;
	}

	Float getSurfaceArea() const {
		Log(EError, "Instance::getSurfaceArea(): not supported!");
		return 0.0f;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->getName() == "ShapeGroup") {
			m_shapeGroup = static_cast<ShapeGroup *>(child);
		} else {
			Shape::addChild(name, child);
		}
	}

	bool rayIntersect(const Ray &_ray, Float mint, 
			Float maxt, Float &t, void *temp) const {
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		Ray ray;
		m_worldToObject(_ray, ray);
		return kdtree->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		Ray ray;
		m_worldToObject(_ray, ray);
		return kdtree->rayIntersect(ray, mint, maxt);
	}

	void fillIntersectionRecord(const Ray &ray, 
		const void *temp, Intersection &its) const {
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		kdtree->fillIntersectionRecord<false>(ray, temp, its);
		its.shFrame.n = m_objectToWorld(its.shFrame.n).normalized();
		its.shFrame.s = m_objectToWorld(its.shFrame.s).normalized();
		its.shFrame.t = m_objectToWorld(its.shFrame.t).normalized();
		its.geoFrame = Frame(m_objectToWorld(its.geoFrame.n).normalized());
		its.wi = its.shFrame.toLocal(-ray.d);
		its.dpdu = m_objectToWorld(its.dpdu);
		its.dpdv = m_objectToWorld(its.dpdv);
	}

	MTS_DECLARE_CLASS()
private:
	ref<ShapeGroup> m_shapeGroup;
	Transform m_objectToWorld, m_worldToObject;
};

MTS_IMPLEMENT_CLASS_S(Instance, false, Shape)
MTS_EXPORT_PLUGIN(Instance, "Instanced geometry");
MTS_NAMESPACE_END
