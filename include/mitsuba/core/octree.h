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

#if !defined(__OCTREE_H)
#define __OCTREE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/octree.h>
#include <mitsuba/core/atomic.h>

MTS_NAMESPACE_BEGIN


/**
 * \brief Implements a lock-free singly linked list using
 * atomic operations.
 * 
 * \ingroup libcore
 */
template <typename T> class LockFreeList {
public:
	struct ListItem {
		T value;
		ListItem *next;

		inline ListItem(const T &value) :
			value(value), next(NULL) { }
	};

	inline LockFreeList() : m_head(NULL) {}

	~LockFreeList() {
		ListItem *cur = m_head;
		while (cur) {
			ListItem *next = cur->next;
			delete cur;
			cur = next;
		}
	}

	inline const ListItem *head() const {
		return m_head;
	}

	void append(const T &value) {
		ListItem *item = new ListItem(value);
		ListItem **cur = &m_head;
		while (!atomicCompareAndExchangePtr<ListItem>(cur, item, NULL))
			cur = &((*cur)->next);
	}
private:
	ListItem *m_head;
};

/**
 * \brief Generic multiple-reference octree.
 *
 * Based on the excellent implementation in PBRT. Modifications are 
 * the addition of a bounding sphere query and support for multithreading.
 *
 * \ingroup libcore
 */
template <typename T> class Octree {
public:
	/**
	 * \brief Create a new octree
	 *
	 * By default, the maximum tree depth is set to 16
	 */

	inline Octree(const BoundingBox3 &bbox, int maxDepth = 16) 
	 : m_bbox(bbox), m_maxDepth(maxDepth) {
	}

	/// Insert an item with the specified cell coverage
	inline void insert(const T &value, const BoundingBox3 &coverage) {
		insert(&m_root, m_bbox, value, coverage,
			coverage.getExtents().squaredNorm(), 0);
	}

	/// Execute operator() of <tt>functor</tt> on all records, which potentially overlap <tt>p</tt>
	template <typename Functor> inline void lookup(const Point &p, Functor &functor) const {
		if (!m_bbox.contains(p))
			return;
		lookup(&m_root, m_bbox, p, functor);
	}

	/// Execute operator() of <tt>functor</tt> on all records, which potentially overlap <tt>bsphere</tt>
	template <typename Functor> inline void searchSphere(const BoundingSphere &sphere, Functor &functor) {
		if (!m_bbox.overlaps(sphere))
			return;
		searchSphere(&m_root, m_bbox, sphere, functor);
	}
	
	inline const BoundingBox3 &getBoundingBox() const { return m_bbox; }
private:
	struct OctreeNode {
	public:
		OctreeNode() {
			for (int i=0; i<8; ++i)
				children[i] = NULL;
		}

		~OctreeNode() {
			for (int i=0; i<8; ++i) {
				if (children[i])
					delete children[i];
			}
		}

		OctreeNode *children[8];
		LockFreeList<T> data;
	};

	/// Return the BoundingBox3 for a child of the specified index
	inline BoundingBox3 childBounds(int child, const BoundingBox3 &nodeBoundingBox, const Point &center) const {
		BoundingBox3 childBoundingBox3;
		childBoundingBox3.min.x = (child & 4) ? center.x : nodeBoundingBox.min.x;
		childBoundingBox3.max.x = (child & 4) ? nodeBoundingBox.max.x : center.x;
		childBoundingBox3.min.y = (child & 2) ? center.y : nodeBoundingBox.min.y;
		childBoundingBox3.max.y = (child & 2) ? nodeBoundingBox.max.y : center.y;
		childBoundingBox3.min.z = (child & 1) ? center.z : nodeBoundingBox.min.z;
		childBoundingBox3.max.z = (child & 1) ? nodeBoundingBox.max.z : center.z;
		return childBoundingBox3;
	}


	void insert(OctreeNode *node, const BoundingBox3 &nodeBoundingBox, const T &value, 
			const BoundingBox3 &coverage, Float diag2, int depth) {
		/* Add the data item to the current octree node if the max. tree
		   depth is reached or the data item's coverage area is smaller
		   than the current node size */
		if (depth == m_maxDepth || 
			(nodeBoundingBox.getExtents().squaredNorm() < diag2)) {
			node->data.append(value);
			return;
		}

		/* Otherwise: test for overlap */
		const Point center = nodeBoundingBox.getCenter();

		/* Otherwise: test for overlap */
		bool x[2] = { coverage.min.x <= center.x, coverage.max.x > center.x };
		bool y[2] = { coverage.min.y <= center.y, coverage.max.y > center.y };
		bool z[2] = { coverage.min.z <= center.z, coverage.max.z > center.z };
		bool over[8] = { x[0] & y[0] & z[0], x[0] & y[0] & z[1],
						 x[0] & y[1] & z[0], x[0] & y[1] & z[1],
						 x[1] & y[0] & z[0], x[1] & y[0] & z[1],
						 x[1] & y[1] & z[0], x[1] & y[1] & z[1] };

		/* Recurse */
		for (int child=0; child<8; ++child) {
			if (!over[child])
				continue;
			if (!node->children[child]) {
				OctreeNode *newNode = new OctreeNode();
				if (!atomicCompareAndExchangePtr<OctreeNode>(&node->children[child], newNode, NULL))
					delete newNode;
			}
			const BoundingBox3 childBoundingBox3(childBounds(child, nodeBoundingBox, center));
			insert(node->children[child], childBoundingBox3,
				value, coverage, diag2, depth+1);
		}
	}

	/// Internal lookup procedure - const version
	template <typename Functor> inline void lookup(const OctreeNode *node, 
			const BoundingBox3 &nodeBoundingBox, const Point &p, Functor &functor) const {
		const Point center = nodeBoundingBox.getCenter();

		const typename LockFreeList<T>::ListItem *item = node->data.head();
		while (item) {
			functor(item->value);
			item = item->next;
		}

		int child = (p.x > center.x ? 4 : 0)
				+ (p.y > center.y ? 2 : 0) 
				+ (p.z > center.z ? 1 : 0);

		OctreeNode *childNode = node->children[child];

		if (childNode) {
			const BoundingBox3 childBoundingBox3(childBounds(child, nodeBoundingBox, center));
			lookup(node->children[child], childBoundingBox3, p, functor);
		}
	}

	template <typename Functor> inline void searchSphere(OctreeNode *node, 
			const BoundingBox3 &nodeBoundingBox, const BoundingSphere &sphere, 
			Functor &functor) {
		const Point center = nodeBoundingBox.getCenter();

		const typename LockFreeList<T>::ListItem *item = node->data.head();
		while (item) {
			functor(item->value);
			item = item->next;
		}

		// Potential for much optimization..
		for (int child=0; child<8; ++child) { 
			if (node->children[child]) {
				const BoundingBox3 childBoundingBox3(childBounds(child, nodeBoundingBox, center));
				if (childBoundingBox3.overlaps(sphere))
					searchSphere(node->children[child], childBoundingBox3, sphere, functor);
			}
		}
	}
private:
	OctreeNode m_root;
	BoundingBox3 m_bbox;
	int m_maxDepth;
};

MTS_NAMESPACE_END

#endif /* __OCTREE_H */
