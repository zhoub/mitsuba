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

#if !defined(__BSPHERE_H)
#define __BSPHERE_H

#include <mitsuba/core/ray.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic n-dimensional bounding sphere data structure
 *
 * This class is parameterized by the underlying point data structure,
 * which permits the use of different scalar types and dimensionalities, e.g.
 * \code
 * BoundingSphere<Point2> doubleBBox(Vector2(0.0f, 1.0f), 2.0f);
 * \endcode
 *
 * \tparam T The underlying point data type (e.g. \c Point2)
 * \sa BoundingSphere
 * \ingroup libcore
 */

template <typename _PointType> struct BoundingSphere {
	typedef _PointType                                   PointType;

	enum {
		Dimension = PointType::RowsAtCompileTime
	};

	typedef typename PointType::Scalar                   Scalar;
	typedef typename Eigen::Matrix<Scalar, Dimension, 1> VectorType;

	PointType center;
	Scalar radius;

	/// Construct a bounding sphere at the origin having radius zero
	inline BoundingSphere() : center(PointType::Zero()), radius(0.0f) { }

	/// Unserialize a bounding sphere from a binary data stream
	inline BoundingSphere(Stream *stream) {
		center = PointType(stream);
		radius = stream->readElement<Scalar>();
	}

	/// Create a bounding sphere from a given center point and radius
	inline BoundingSphere(const PointType &center, Scalar radius)
		: center(center), radius(radius) { }

	/// Return whether this bounding sphere has a radius of zero or less.
	inline bool isEmpty() const {
		return radius <= 0.0f;
	}

	/// Return the bounding sphere's center 
	inline const PointType &getCenter() const { return center; }

	/// Return the bounding sphere's radius
	inline Scalar getRadius() const { return radius; }

	/// Expand the bounding sphere radius to contain another point.
	inline void expandBy(const PointType p) {
		radius = std::max(radius, (p-center).norm());
	}

	/**
	 * \brief Check whether the specified point lies \a inside or \a on the sphere
	 *
	 * \param p The point to be tested
	 *
	 * \param strict Set this parameter to \c true if the bounding
	 *               sphere boundary should be excluded in the test
	 */
	inline bool contains(const PointType p, bool strict = false) const {
		if (strict)
			return (p - center).norm() < radius;
		else
			return (p - center).norm() <= radius;
	}

	/// Test for equality against another bounding sphere
	inline bool operator==(const BoundingSphere &boundingSphere) const {
		return center == boundingSphere.center && radius == boundingSphere.radius;
	}

	/// Test for inequality against another bounding sphere
	inline bool operator!=(const BoundingSphere &boundingSphere) const {
		return center != boundingSphere.center || radius != boundingSphere.radius;
	}

	/// Serialize this bounding sphere to a binary data stream
	inline void serialize(Stream *stream) const {
		center.serialize(stream);
		stream->writeElement(radius);
	}

	/// Return a string representation of the bounding sphere
	inline std::string toString() const {
		std::ostringstream oss;
		oss << "BoundingSphere[center = " << center.toString()
			<< ", radius = " << radius << "]";
		return oss.str();
	}
};

/**
 * \brief Bounding sphere data structure in three dimensions
 * 
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE BoundingSphere3 : public BoundingSphere<Point> {
	/// Construct a bounding sphere at the origin having radius zero
	inline BoundingSphere3() : BoundingSphere<Point>() { }

	/// Unserialize a bounding sphere from a binary data stream
	inline BoundingSphere3(Stream *stream) : BoundingSphere<Point>(stream) { }

	/// Create a bounding sphere from a given center point and radius
	inline BoundingSphere3(const Point &center, Float radius)
		: BoundingSphere<Point>(center, radius) { }

	/// Construct from a BoundingSphere<Point>
	inline BoundingSphere3(const BoundingSphere<Point> &bsphere) 
		: BoundingSphere<Point>(bsphere) { }

	/**
	 * \brief Calculate the intersection points with the given ray
	 *
	 * Also returns intersection points along the negative ray direction.
	 *
	 * \return \c true if the ray intersects the bounding sphere
	 */
	inline bool rayIntersect(const Ray &ray, Scalar &nearHit, Scalar &farHit) const {
		const Vector o = ray.o - center;
		const Float A = ray.d.squaredNorm();
		const Float B = 2 * ray.d.dot(o);
		const Float C = o.squaredNorm() - radius*radius;

		if (!solveQuadratic(A, B, C, nearHit, farHit))
			return false;

		return true;
	}
};


MTS_NAMESPACE_END

#endif /* __BSPHERE_H */
