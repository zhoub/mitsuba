/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#pragma once
#if !defined(__MITSUBA_RENDER_VKDTREE2_H_)
#define __MITSUBA_RENDER_VKDTREE2_H_

#include <mitsuba/core/aabb.h>
#include <mitsuba/render/gkdtree.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Implements a 2D "volume" heuristic for use
 * by the \ref GenericKDTree construction algorithm.
 * \ingroup librender
 */
class VolumeHeuristic2 {
public:
	/**
	 * \brief Initialize the surface area heuristic with the bounds of
	 * a parent node
	 *
	 * Precomputes some information so that traversal probabilities
	 * of potential split planes can be evaluated efficiently
	 */
	inline VolumeHeuristic2(const AABB2 &aabb) {
		m_extents = aabb.getExtents();
		m_normalization = 1.0f / (m_extents.x * m_extents.y);
	}

	/**
	 * Given a split on axis \a axis that produces children having extents
	 * \a leftWidth and \a rightWidth along \a axis, compute the probability
	 * of traversing the left and right child during a typical query
	 * operation.
	 */
	inline std::pair<Float, Float> operator()(int axis, Float leftWidth, Float rightWidth) const {
		return std::pair<Float, Float>(
			(m_extents[1-axis] * leftWidth) * m_normalization,
			(m_extents[1-axis] * rightWidth) * m_normalization);
	}

	/**
	 * Compute the underlying quantity used by the tree construction
	 * heuristic. This is used to compute the final cost of a kd-tree.
	 */
	inline static Float getQuantity(const AABB2 &aabb) {
		return aabb.getVolume();
	}
private:
	Float m_normalization;
	Vector2 m_extents;
};

/**
 * \brief Specializes \ref GenericKDTree to a two-dimensional
 * tree to be used for parameterizing shapes
 *
 * One additional function call must be implemented by subclasses:
 * \code
 * /// Check whether a primitive contains the specified point
 * bool contains(IndexType idx, const Point2 &p, Point2 &uv);
 * \endcode
 *
 * \author Wenzel Jakob
 * \ingroup librender
 */
template <typename Derived>
	class VKDTree2D : public GenericKDTree<AABB2, VolumeHeuristic2, Derived> {
public:
	typedef GenericKDTree<AABB2, VolumeHeuristic2, Derived>      Parent;
	typedef typename KDTreeBase<AABB2>::SizeType                 SizeType;
	typedef typename KDTreeBase<AABB2>::IndexType                IndexType;
	typedef typename KDTreeBase<AABB2>::KDNode                   KDNode;

	using Parent::m_nodes;
	using Parent::m_aabb;
	using Parent::m_indices;

protected:
	void buildInternal() {
		SizeType primCount = cast()->getPrimitiveCount();
		KDLog(EInfo, "Constructing a 2D kd-tree (%i primitives) ..", primCount);
		GenericKDTree<AABB2, VolumeHeuristic2, Derived>::buildInternal();
	}

	/// Cast to the derived class
	inline Derived *cast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *cast() const {
		return static_cast<const Derived *>(this);
	}

	FINLINE bool find(const Point2 &p, uint32_t &index, Point2 &uv) const {
		if (!m_aabb.contains(p))
			return false;

		const KDNode * __restrict node = m_nodes;
		while (!node->isLeaf()) {
			const KDNode * __restrict left = node->getLeft();
			node = left + (p[node->getAxis()] < (Float) node->getSplit() ? 0 : 1);
		}

		for (uint32_t entry=node->getPrimStart(),
				last = node->getPrimEnd(); entry != last; entry++) {
			const IndexType primIdx = m_indices[entry];

			if (cast()->contains(primIdx, p, uv)) {
				index = primIdx;
				return true;
			}
		}
		return false;
	}
};

/**
 * \brief KD-tree class used to parameterize 3D triangle shapes
 * using a 2D domain
 */
class UVKDTree : public VKDTree2D<UVKDTree> {
	friend class GenericKDTree<AABB2, VolumeHeuristic2, UVKDTree>;
	friend class VKDTree2D<UVKDTree>;
public:
	using VKDTree2D<UVKDTree>::IndexType;
	using VKDTree2D<UVKDTree>::SizeType;
	using VKDTree2D<UVKDTree>::find;

	UVKDTree(SizeType triangleCount, const Triangle *triangles, const Point2 *uv)
		 : m_triangleCount(triangleCount), m_triangles(triangles), m_uv(uv) {
		buildInternal();
	}

	AABB2 getAABB(IndexType index) const {
		const Triangle &tri = m_triangles[index];
		AABB2 aabb(m_uv[tri.idx[0]]);

		for (int i=1; i<3; ++i)
			aabb.expandBy(m_uv[tri.idx[i]]);

		return aabb;
	}

	/// Compute the clipped AABB of a segment (only used during tree construction)
	AABB2 getClippedAABB(IndexType index, const AABB2 &box) const {
		AABB2 aabb(getAABB(index)); /* Could do a better job here */
		aabb.clip(box);
		return aabb;
	}

	bool contains(IndexType index, const Point2 &p, Point2 &bary) const {
		const Triangle &tri = m_triangles[index];
		const Point2 &p0 = m_uv[tri.idx[0]];
		Vector2 rel = p - p0;
		Vector2 d1 = m_uv[tri.idx[1]] - p0;
		Vector2 d2 = m_uv[tri.idx[2]] - p0;

		Float det = d1.y*d2.x-d1.x*d2.y;
		if (det == 0)
			return false;

		Float invDet = (Float) 1 / det;
		Float b0 = (-d2.y*rel.x + d2.x*rel.y) * invDet;
		Float b1 = ( d1.y*rel.x - d1.x*rel.y) * invDet;

		if (!(b0 >= 0 && b0 <= 1 && b1 >= 0 && b1 <= 1 && b0 + b1 <= 1))
			return false;

		bary = Point2(b0, b1);
		return true;
	}

	/// Return the total number of segments
	inline SizeType getPrimitiveCount() const {
		return m_triangleCount;
	}

	MTS_DECLARE_CLASS()
protected:
	SizeType m_triangleCount;
	const Triangle *m_triangles;
	const Point2 *m_uv;
};


MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_VKDTREE2_H_ */
