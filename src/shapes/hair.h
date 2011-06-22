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

#if !defined(__HAIR_H)
#define __HAIR_H

#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN
	
class HairKDTree;

/**
 * \brief Intersection shape structure for cylindrical hair
 * segments with miter joints. This class expects an ASCII file containing
 * a list of hairs made from segments. Each line should contain an X,
 * Y and Z coordinate separated by a space. An empty line indicates
 * the start of a new hair.
 */
class HairShape : public Shape {
public:
	/// Construct a new HairShape instance given a properties object
	HairShape(const Properties &props);
	
	/// Unserialize from a binary data stream
	HairShape(Stream *stream, InstanceManager *manager);

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	// =============================================================
	//! @{ \name Access the internal vertex data
	// =============================================================
	
	/// Return the list of vertices underlying the hair shape
	//const std::vector<Point> &getVertices() const;

	/**
	 * Return a boolean list specifying whether a vertex 
	 * marks the beginning of a new fiber
	 */
	//const std::vector<bool> &getStartFiber() const;

	// Get the radius of the fibers
	Float getRadius() const { return (this->*ptr_getRadius)(); }

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Geometric information about hair segments
	// =============================================================

	// Extract information about hit location that is stored in the
	// intersection record at intersection time.  Output is
	// the index of the segment that was hit and the (u,v,s)
	// segment coordinates of the hit point.
	void convertLocationData(const Intersection &its, int &iSeg,
			Point &segPosn) const {
		(this->*ptr_convertLocationData)(its, iSeg, segPosn);
	}

	// Convert a point in (u,v,s) segment coordinates to world space.
	Point segmentToGlobal(const Intersection &its, int iSeg,
			const Point &segPosn) const {
		return (this->*ptr_segmentToGlobal)(its, iSeg, segPosn);
	}

	// The first and second vertices of a given segment
	Point getFirstVertex(int iSeg) const {
		return (this->*ptr_getFirstVertex)(iSeg);
	}
	Point getSecondVertex(int iSeg) const {
		return (this->*ptr_getSecondVertex)(iSeg);
	}

	// The tangent vector of a given segment
	Vector getSegmentTangent(int iSeg) const {
		return (this->*ptr_getSegmentTangent)(iSeg);
	}

	// The miter normals of a given segment
	Vector getFirstMiterNormal(int iSeg) const {
		return (this->*ptr_getFirstMiterNormal)(iSeg);
	}
	Vector getSecondMiterNormal(int iSeg) const {
		return (this->*ptr_getSecondMiterNormal)(iSeg);
	}

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Implementation of the \ref Shape interface
	// =============================================================

	bool rayIntersect(const Ray &ray, Float mint, 
			Float maxt, Float &t, void *temp) const;
	
	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const;

	void fillIntersectionRecord(const Ray &ray, 
		const void *temp, Intersection &its) const;

	ref<TriMesh> createTriMesh();

	const KDTreeBase<AABB> *getKDTree() const;

	AABB getAABB() const;

	Float getSurfaceArea() const;
	
	//! @}
	// =============================================================

	/// Return a human-readable representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
private:
	ref<HairKDTree> m_kdtree;

	Float (HairShape::*ptr_getRadius)() const;
	void (HairShape::*ptr_convertLocationData)(const Intersection &, int &,
			Point &) const;
	Point (HairShape::*ptr_segmentToGlobal)(const Intersection &, int,
			const Point &) const;
	Point (HairShape::*ptr_getFirstVertex)(int) const;
	Point (HairShape::*ptr_getSecondVertex)(int) const;
	Vector (HairShape::*ptr_getSegmentTangent)(int) const;
	Vector (HairShape::*ptr_getFirstMiterNormal)(int) const;
	Vector (HairShape::*ptr_getSecondMiterNormal)(int) const;

	Float impl_getRadius() const;
	void impl_convertLocationData(const Intersection &its, int &iSeg,
			Point &segPosn) const;
	Point impl_segmentToGlobal(const Intersection &its, int iSeg,
			const Point &segPosn) const;
	Point impl_getFirstVertex(int iSeg) const;
	Point impl_getSecondVertex(int iSeg) const;
	Vector impl_getSegmentTangent(int iSeg) const;
	Vector impl_getFirstMiterNormal(int iSeg) const;
	Vector impl_getSecondMiterNormal(int iSeg) const;

	void initImplementations();

};

MTS_NAMESPACE_END

#endif /* __HAIR_H */
