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

#include <mitsuba/render/track.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include "shapegroup.h"

MTS_NAMESPACE_BEGIN

class AnimatedInstance : public Shape {
public:
	AnimatedInstance(const Properties &props) : Shape(props) {
		FileResolver *fResolver = Thread::getThread()->getFileResolver();
		fs::path path = fResolver->resolve(props.getString("filename"));
		m_name = path.filename();

		Log(EInfo, "Loading animation track from \"%s\"", m_name.c_str());
		ref<FileStream> fs = new FileStream(path, FileStream::EReadOnly);
		m_transform = new AnimatedTransform(fs);
	}

	AnimatedInstance(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		m_shapeGroup = static_cast<ShapeGroup *>(manager->getInstance(stream));
		m_transform = new AnimatedTransform(stream);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		manager->serialize(stream, m_shapeGroup.get());
		m_transform->serialize(stream);
	}

	void configure() {
		if (!m_shapeGroup)
			Log(EError, "A reference to a 'shapegroup' must be specified!");
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		const BoundingBox3 &bbox = kdtree->getBoundingBox();
		Float minT, maxT;
		m_transform->computeTimeBounds(minT, maxT);

		/* Compute approximate bounds */
		int nSteps = 100;
		Float step = (maxT-minT) / (nSteps-1);
		Transform objectToWorld;

		for (int i=0; i<nSteps; ++i) {
			m_transform->eval(minT + step * i, objectToWorld);
			for (int j=0; j<8; ++j)
				m_bbox.expandBy(objectToWorld(bbox.getCorner(j)));
		}
	}

	BoundingBox3 getBoundingBox() const {
		return m_bbox;
	}

	std::string getName() const {
		return m_name;
	}

	Float getSurfaceArea() const {
		Log(EError, "AnimatedInstance::getSurfaceArea(): not supported!");
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
		Transform objectToWorld, worldToObject;
		m_transform->eval(_ray.time, objectToWorld);
		worldToObject = objectToWorld.inverse();
		worldToObject(_ray, ray);
		return kdtree->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		Ray ray;
		Transform objectToWorld, worldToObject;
		m_transform->eval(_ray.time, objectToWorld);
		worldToObject = objectToWorld.inverse();
		worldToObject(_ray, ray);
		return kdtree->rayIntersect(ray, mint, maxt);
	}

	void fillIntersectionRecord(const Ray &ray, 
		const void *temp, Intersection &its) const {
		const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
		Transform objectToWorld;
		m_transform->eval(ray.time, objectToWorld);
		kdtree->fillIntersectionRecord<false>(ray, temp, its);
		its.shFrame.n = (objectToWorld(its.shFrame.n)).normalized();
		its.shFrame.s = (objectToWorld(its.shFrame.s)).normalized();
		its.shFrame.t = (objectToWorld(its.shFrame.t)).normalized();
		its.geoFrame = Frame(objectToWorld(its.geoFrame.n).normalized());
		its.wi = its.shFrame.toLocal(-ray.d);
		its.dpdu = objectToWorld(its.dpdu);
		its.dpdv = objectToWorld(its.dpdv);
	}

	MTS_DECLARE_CLASS()
private:
	ref<ShapeGroup> m_shapeGroup;
	ref<AnimatedTransform> m_transform;
	BoundingBox3 m_bbox;
	std::string m_name;
};

MTS_IMPLEMENT_CLASS_S(AnimatedInstance, false, Shape)
MTS_EXPORT_PLUGIN(AnimatedInstance, "Animated instanced geometry");
MTS_NAMESPACE_END
