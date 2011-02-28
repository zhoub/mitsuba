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

#if !defined(__BOUNDING_BOX_H)
#define __BOUNDING_BOX_H

#include <mitsuba/core/bsphere.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic n-dimensional bounding box data structure
 *
 * Maintains a minimum and maximum position along each dimension and provides
 * various convenience functions for querying and modifying them.
 *
 * This class is parameterized by the underlying point data structure,
 * which permits the use of different scalar types and dimensionalities, e.g.
 * \code
 * BoundingBox<Vector3i> integerBBox(Point3i(0, 1, 3), Point3i(4, 5, 6));
 * BoundingBox<Vector2d> doubleBBox(Point2d(0.0, 1.0), Point2d(4.0, 5.0));
 * \endcode
 *
 * \tparam T The underlying point data type (e.g. \c Point2d)
 * \sa BoundingSphere
 * \ingroup libcore
 */
template <typename _PointType> struct BoundingBox {
	typedef _PointType                                   PointType;

	enum {
		Dimension = PointType::RowsAtCompileTime
	};

	typedef typename PointType::Scalar                   Scalar;
	typedef typename Eigen::Matrix<Scalar, Dimension, 1> VectorType;

	/** 
	 * \brief Create a new invalid bounding box
	 * 
	 * Initializes the components of the minimum 
	 * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	 * respectively.
	 */
	inline BoundingBox() {
		reset();
	}

	/// Unserialize a bounding box from a binary data stream
	inline BoundingBox(Stream *stream) {
		min = PointType(stream);
		max = PointType(stream);
	}

	/// Create a collapsed bounding box from a single point
	inline BoundingBox(const PointType &p) 
		: min(p), max(p) { }

	/// Create a bounding box from two positions
	inline BoundingBox(const PointType &min, const PointType &max)
		: min(min), max(max) {
#if defined(MTS_DEBUG)
		SAssert(isValid());
#endif
	}

	/// Test for equality against another bounding box
	inline bool operator==(const BoundingBox &bbox) const {
		return min == bbox.min && max == bbox.max;
	}

	/// Test for inequality against another bounding box
	inline bool operator!=(const BoundingBox &bbox) const {
		return min != bbox.min || max != bbox.max;
	}

	/// Calculate the n-dimensional volume of the bounding box
	inline Scalar getVolume() const {
		return (max - min).prod();
	}

	/// Calculate the n-1 dimensional volume of the boundary
	inline Float getSurfaceArea() const {
		VectorType d = max - min;
		Float result = 0.0f;
		for (int i=0; i<Dimension; ++i) {
			Float term = 1.0f;
			for (int j=0; j<Dimension; ++j) {
				if (i == j)
					continue;
				term *= d[j];
			}
			result += term;
		}
		return 2.0f * result;
	}

	/// Return the center point
	inline PointType getCenter() const {
		return (max + min) * (Scalar) 0.5f;
	}

	/**
	 * \brief Check whether a point lies \a on or \a inside the bounding box
	 *
	 * \param p The point to be tested
	 *
	 * \param strict Set this parameter to \c true if the bounding
	 *               box boundary should be excluded in the test
	 */
	inline bool contains(const PointType &p, bool strict = false) const {
		if (strict) {
			return (p.array() > min.array()).all() 
				&& (p.array() < max.array()).all();
		} else {
			return (p.array() >= min.array()).all() 
				&& (p.array() <= max.array()).all();
		}
	}

	/**
	 * \brief Check whether a specified bounding box lies \a on or \a within 
	 * the current bounding box
	 *
	 * Note that by definition, an 'invalid' bounding box (where min=\f$\infty\f$
	 * and max=\f$-\infty\f$) does not cover any space. Hence, this method will always 
	 * return \a true when given such an argument.
	 *
	 * \param strict Set this parameter to \c true if the bounding
	 *               box boundary should be excluded in the test
	 */
	inline bool contains(const BoundingBox &bbox, bool strict = false) const {
		if (strict) {
			return (bbox.min.array() > min.array()).all() 
				&& (bbox.max.array() < max.array()).all();
		} else {
			return (bbox.min.array() >= min.array()).all() 
				&& (bbox.max.array() <= max.array()).all();
		}
	}

	/**
	 * \brief Check two axis-aligned bounding boxes for possible overlap.
	 *
	 * \param strict Set this parameter to \c true if the bounding
	 *               box boundary should be excluded in the test
	 *
	 * \return \c true If overlap was detected.
	 */
	inline bool overlaps(const BoundingBox &bbox, bool strict = false) const {
		if (strict) {
			return (bbox.min.array() < max.array()).all() 
				&& (bbox.max.array() > min.array()).all();
		} else {
			return (bbox.min.array() <= max.array()).all() 
				&& (bbox.max.array() >= min.array()).all();
		}
	}

	/**
	 * \brief Calculate the smallest squared distance between
	 * the axis-aligned bounding box and the point \c p.
	 */
	inline Scalar squaredDistanceTo(const PointType &p) const {
		Scalar result = 0;

		for (int i=0; i<Dimension; ++i) {
			Scalar value = 0;
			if (p[i] < min[i])
				value = min[i] - p[i];
			else if (p[i] > max[i])
				value = p[i] - max[i];
			result += value*value;
		}

		return result;
	}

	/**
	 * \brief Calculate the smallest distance between
	 * the axis-aligned bounding box and the point \c p.
	 */
	inline Scalar distanceTo(const PointType &p) const {
		return std::sqrt(squaredDistanceTo(p));
	}

	/**
	 * \brief Calculate the smallest square distance between
	 * the axis-aligned bounding box and \c bbox.
	 */
	inline Scalar squaredDistanceTo(const BoundingBox &bbox) const {
		Scalar result = 0;

		for (int i=0; i<Dimension; ++i) {
			Scalar value = 0;
			if (bbox.max[i] < min[i])
				value = min[i] - bbox.max[i];
			else if (bbox.min[i] > max[i])
				value = bbox.min[i] - max[i];
			result += value*value;
		}

		return result;
	}

	/**
	 * \brief Calculate the smallest distance between
	 * the axis-aligned bounding box and \c bbox.
	 */
	inline Scalar distanceTo(const BoundingBox &bbox) const {
		return std::sqrt(squaredDistanceTo(bbox));
	}

	/**
	 * \brief Check whether this is a valid bounding box
	 *
	 * A bounding box \c bbox is valid when
	 * \code
	 * bbox.min[dim] <= bbox.max[dim]
	 * \endcode
	 * holds along each dimension \c dim.
	 */
	inline bool isValid() const {
		return (max.array() >= min.array()).all();
	}

	/// Check whether this bounding box has collapsed to a single point
	inline bool isPoint() const {
		return (max.array() == min.array()).all();
	}

	/// Check whether this bounding box has any associated volume
	inline bool hasVolume() const {
		return (max.array() > min.array()).all();
	}

	/// Return the dimension index with the largest associated side length
	inline int getMajorAxis() const {
		VectorType d = max - min;
		int largest = 0;
		for (int i=1; i<Dimension; ++i)
			if (d[i] > d[largest])
				largest = i;
		return largest;
	}

	/// Return the dimension index with the shortest associated side length
	inline int getMinorAxis() const {
		VectorType d = max - min;
		int shortest = 0;
		for (int i=1; i<Dimension; ++i)
			if (d[i] < d[shortest])
				shortest = i;
		return shortest;
	}

	/**
	 * \brief Calculate the bounding box extents
	 * \return max-min
	 */
	inline VectorType getExtents() const {
		return max - min;
	}

	/// Create a bounding sphere, which contains the axis-aligned box
	inline BoundingSphere<PointType> getBoundingSphere() const {
		PointType center = getCenter();
		return BoundingSphere<PointType>(center, (center - max).norm());
	}

	/// Serialize this bounding box to a binary data stream
	inline void serialize(Stream *stream) const {
		min.serialize(stream);
		max.serialize(stream);
	}

	/// Clip to another bounding box
	inline void clip(const BoundingBox &bbox) {
		min = min.cwiseMax(bbox.min);
		max = max.cwiseMin(bbox.max);
	}

	/** 
	 * \brief Mark the bounding box as invalid.
	 * 
	 * This operation sets the components of the minimum 
	 * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	 * respectively.
	 */
	inline void reset() {
		min.setConstant( std::numeric_limits<Scalar>::infinity());
		max.setConstant(-std::numeric_limits<Scalar>::infinity());
	}

	/// Expand the bounding box to contain another point
	inline void expandBy(const PointType &p) {
		min = min.cwiseMin(p);
		max = max.cwiseMax(p);
	}

	/// Expand the bounding box to contain another bounding box
	inline void expandBy(const BoundingBox &bbox) {
		min = min.cwiseMin(bbox.min);
		max = max.cwiseMax(bbox.max);
	}

	/// Return a string representation of the bounding box
	std::string toString() const {
		std::ostringstream oss;
		oss << "BoundingBox3[";
		if (!isValid()) {
			oss << "invalid";
		} else {
			oss << "min=" << min.toString()
				<< ", max=" << max.toString();
		}
		oss	<< "]";
		return oss.str();
	}

	PointType min; ///< Component-wise minimum 
	PointType max; ///< Component-wise maximum 
};


/**
 * \brief Axis-aligned bounding box data structure in three dimensions
 * 
 * Maintains a component-wise minimum and maximum position and provides
 * various convenience functions to query or change them.
 *
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE BoundingBox3 : public BoundingBox<Point> {
public:
	using BoundingBox<Point>::overlaps;

	/** 
	 * \brief Create a new invalid bounding box
	 * 
	 * Initializes the components of the minimum 
	 * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	 * respectively.
	 */
	inline BoundingBox3() : BoundingBox<Point>() { }

	/// Unserialize a bounding box from a binary data stream
	inline BoundingBox3(Stream *stream) : BoundingBox<Point>(stream) { }

	/// Create a collapsed bounding box from a single point
	inline BoundingBox3(const Point &p) : BoundingBox<Point>(p) { }

	/// Create a bounding box from two positions
	inline BoundingBox3(const PointType &min, const PointType &max) 
		: BoundingBox<Point>(min, max) {
	}

	/// Construct from a BoundingBox<Point>
	inline BoundingBox3(const BoundingBox<Point> &bbox) 
		: BoundingBox<Point>(bbox) { }

	/// Calculate the surface area of the bounding box
	inline Float getSurfaceArea() const {
		Vector d = max - min;
		return (Float) 2.0f * (d[0]*d[1] + d[0]*d[2] + d[1]*d[2]);
	}

	/**
	 * \brief Return the position of a bounding box corner
	 * \param corner Requested corner index (0..7)
	 */
	inline Point getCorner(int corner) const {
		return Point(
				corner & 1 ? max[0] : min[0],
				corner & 2 ? max[1] : min[1],
				corner & 4 ? max[2] : min[2]);
	}

	/**
	 * \brief Bounding sphere-box overlap test
	 *
	 * Implements the technique proposed by Jim Arvo in
	 * "A simple method for box-sphere intersection testing"
	 * (Graphics Gems, 1990)
	 */
	bool overlaps(const BoundingSphere3 &sphere) const {
		Float distance = 0;
		for (int i=0; i<3; ++i) {
			if (sphere.center[i] < min[i]) {
				Float d = sphere.center[i]-min[i];
				distance += d*d;
			} else if (sphere.center[i] > max[i]) {
				Float d = sphere.center[i]-max[i];
				distance += d*d;
			}
		}
		return distance < sphere.radius*sphere.radius;
	}

	/** \brief Calculate the near and far ray-box intersection
	 * points (if they exist).
	 *
	 * Also returns intersections along the negative ray direction.
	 */
	FINLINE bool rayIntersect(const Ray &ray, Float &nearT, Float &farT) const {
		nearT = -std::numeric_limits<Float>::infinity();
		farT  = std::numeric_limits<Float>::infinity();

		/* For each pair of bounding box planes */
		for (int i=0; i<3; i++) {
			const Float direction = ray.d[i];
			const Float origin = ray.o[i];
			const Float minVal = min[i], maxVal = max[i];

			if (direction == 0) {
				/* The ray is parallel to the planes */
				if (origin < minVal || origin > maxVal)
					return false;
			} else {
				/* Calculate intersection distances */
				Float t1 = (minVal - origin) * ray.dRcp[i];
				Float t2 = (maxVal - origin) * ray.dRcp[i];

				if (t1 > t2) 
					std::swap(t1, t2);

				nearT = std::max(nearT, t1);
				farT = std::min(farT, t2);

				if (nearT > farT)
					return false;
			}
		}
		return true;
	}

#ifdef MTS_SSE
	/**
	 * \brief Intersect against a packet of four rays. 
	 * \return \a false if none of the rays intersect.
	 */
	FINLINE bool rayIntersectPacket(const RayPacket4 &ray, RayInterval4 &interval) const;
#endif
};

MTS_NAMESPACE_END

#ifdef MTS_SSE
#include <mitsuba/core/bbox_sse.h>
#endif

#endif /* __BOUNDING_BOX_H */
