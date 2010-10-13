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

#if !defined(__KD_TREE_H)
#define __KD_TREE_H

#include <mitsuba/render/shape.h>
#include <mitsuba/render/gkdtree.h>
#include <mitsuba/render/triaccel.h>

#define MTS_KD_CONSERVE_MEMORY 1

MTS_NAMESPACE_BEGIN

/**
 * \brief SAH KD-tree acceleration data structure for fast ray-triangle 
 * intersections.
 *
 * Implements the construction algorithm for 'perfect split' trees as outlined 
 * in the paper "On Bulding fast kd-Trees for Ray Tracing, and on doing that in
 * O(N log N)" by Ingo Wald and Vlastimil Havran. Non-triangle shapes are 
 * supported, but most optimizations here target large triangle meshes.
 * For more details regarding the construction algorithm, please refer to
 * the class \ref GenericKDTree.
 *
 * This class offers a choice of two different triangle intersection algorithms:
 * By default, intersections are computed using the "TriAccel" projection with 
 * pre-computation method from Ingo Wald's PhD thesis "Realtime Ray Tracing 
 * and Interactive Global Illumination". This adds an overhead of 48 bytes per
 * triangle.
 *
 * When compiled with MTS_KD_CONSERVE_MEMORY, the Moeller-Trumbore intersection 
 * test is used instead, which doesn't need any extra storage. However, it also
 * tends to be quite a bit slower.
 *
 * \sa GenericKDTree
 */
class KDTree : public GenericKDTree<KDTree> {
	friend class GenericKDTree<KDTree>;
public:
	/// Create an empty kd-tree
	KDTree();

	/// Add a shape to the kd-tree
	void addShape(const Shape *shape);

	/// Return the list of stored shapes
	inline const std::vector<const Shape *> &getShapes() const { return m_shapes; }

	/**
	 * \brief Return an axis-aligned bounding box containing all primitives
	 */
	inline const AABB &getAABB() const { return m_aabb; }

	/// Build the kd-tree (needs to be called before tracing any rays)
	void build();

	/**
	 * \brief Intersect a ray against all primitives stored in the kd-tree
	 */
	bool rayIntersect(const Ray &ray, Intersection &its) const;

	/**
	 * \brief Test a ray for intersection against all primitives stored in the kd-tree
	 */
	bool rayIntersect(const Ray &ray) const;

	MTS_DECLARE_CLASS()
protected:
	/**
	 * \brief Return the shape index corresponding to a primitive index
	 * seen by the generic kd-tree implementation. When this is a triangle
	 * mesh, the \a idx parameter is updated to the triangle index within
	 * the mesh.
	 */
	index_type findShape(index_type &idx) const {
		std::vector<index_type>::const_iterator it = std::lower_bound(
				m_shapeMap.begin(), m_shapeMap.end(), idx+1) - 1;
		idx -= *it;
		return (index_type) (it - m_shapeMap.begin());
	}

 	/// Return the axis-aligned bounding box of a certain primitive
	FINLINE AABB getAABB(index_type idx) const {
		index_type shapeIdx = findShape(idx);
		const Shape *shape = m_shapes[shapeIdx];
		if (m_triangleFlag[shapeIdx]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			return mesh->getTriangles()[idx].getAABB(mesh->getVertexBuffer());
		} else {
			return shape->getAABB();
		}
	}

 	/// Return the AABB of a primitive when clipped to another AABB
	FINLINE AABB getClippedAABB(index_type idx, const AABB &aabb) const {
		index_type shapeIdx = findShape(idx);
		const Shape *shape = m_shapes[shapeIdx];
		if (m_triangleFlag[shapeIdx]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			return mesh->getTriangles()[idx].getClippedAABB(mesh->getVertexBuffer(), aabb);
		} else {
			return shape->getAABB();
		}
	}

	FINLINE size_type getPrimitiveCount() const {
		return m_shapeMap[m_shapeMap.size()-1];
	}
 
	/// Temporarily holds some intersection information
	struct IntersectionCache {
		size_type shapeIndex;
		size_type index;
		Float u, v;
	};

	/**
	 * Check whether a primitive is intersected by the given ray. Some
	 * temporary space is supplied to store data that can later
	 * be used to create a detailed intersection record.
	 */
	FINLINE EIntersectionResult intersect(const Ray &ray, index_type idx, Float mint, 
		Float maxt, Float &t, void *temp) const {
		Float tempU, tempV, tempT;
	
#if defined(MTS_KD_CONSERVE_MEMORY)
		index_type shapeIdx = findShape(idx);
		if (EXPECT_TAKEN(m_triangleFlag[shapeIdx])) {
			const TriMesh *mesh = static_cast<const TriMesh *>(m_shapes[shapeIdx]);
			const Triangle &tri = mesh->getTriangles()[idx];
			if (tri.rayIntersect(mesh->getVertexBuffer(), ray, tempU, tempV, tempT)) {
				if (tempT < mint || tempT > maxt)
					return ENo;
				IntersectionCache *cache = static_cast<IntersectionCache *>(temp);
				t = tempT;
				cache->shapeIndex = shapeIdx;
				cache->index = idx;
				cache->u = tempU;
				cache->v = tempV;
				return EYes;
			}
		} else {
			cout << "Encountered a non-triangle shape!" << endl;
		}
#else
		if (EXPECT_TAKEN(m_triAccel[idx].k != KNoTriangleFlag)) {
			const TriAccel &ta = m_triAccel[idx];
			if (ta.rayIntersect(ray, mint, maxt, tempU, tempV, tempT)) {
				IntersectionCache *cache = static_cast<IntersectionCache *>(temp);
				t = tempT;
				cache->shapeIndex = ta.shapeIndex;
				cache->index = ta.index;
				cache->u = tempU;
				cache->v = tempV;
				return EYes;
			}
		} else {
			cout << "Encountered a non-triangle shape!" << endl;
			//int shape = m_triAccel[idx].shapeIndex;
		}
#endif
		return ENo;
	}

	/// Virtual destructor
	virtual ~KDTree();
private:
	std::vector<const Shape *> m_shapes;
	std::vector<bool> m_triangleFlag;
	std::vector<index_type> m_shapeMap;
#if !defined(MTS_KD_CONSERVE_MEMORY)
	TriAccel *m_triAccel;
#endif
};

MTS_NAMESPACE_END

#endif /* __KD_TREE_H */
