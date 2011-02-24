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

#if !defined(__VECTOR_H)
#define __VECTOR_H

namespace mitsuba { class Stream; };

#define EIGEN_MATRIX_PLUGIN "mitsuba/core/eigen_matrix.inl"
#define EIGEN_ARRAY_PLUGIN "mitsuba/core/eigen_array.inl"

#include <Eigen/Dense>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic N-dimensional point data structure based on Eigen::Matrix.
 *
 * The sole reason for the existence of this class is that points must be
 * treated differently from vectors when undergoing a transformation in homogeneous
 * coordinates. See \ref Transform for more details.
 */
template <typename T, int N> class TPoint : public Eigen::Matrix<T, N, 1> {
public:
	typedef Eigen::Matrix<T, N, 1> Base;

	inline TPoint() : Base() { }
	inline TPoint(T x, T y) : Base(x, y) { }
	inline TPoint(T x, T y, T z) : Base(x, y, z) { }
	inline TPoint(T x, T y, T z, T w) : Base(x, y, z, w) { }
	inline TPoint(Stream *stream) : Base(stream) { }
	template<typename Derived> inline TPoint(const Eigen::MatrixBase<Derived>& p) 
		: Base(p) { }

    template<typename Derived> TPoint &operator=
			(const Eigen::MatrixBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
    }
};

/**
 * \brief 3D normal vector data structure based on Eigen::Matrix.
 *
 * The sole reason for the existence of this class is that normal vectors must be
 * treated differently from points and vectors when undergoing a transformation 
 * in homogeneous coordinates. See \ref Transform for more details.
 */
class Normal : public Eigen::Matrix<Float, 3, 1> {
public:
	typedef Eigen::Matrix<Float, 3, 1> Base;

	inline Normal() : Base() { }
	inline Normal(Float x, Float y, Float z) : Base(x, y, z) { }
	inline Normal(Stream *stream) : Base(stream) { }
	template<typename Derived> inline Normal(const Eigen::MatrixBase<Derived>& n) 
		: Base(n) { }

    template<typename OtherDerived> Normal &operator=
			(const Eigen::MatrixBase<OtherDerived>& other) {
		this->Base::operator=(other);
		return *this;
    }
};


typedef Eigen::Matrix<Float,  3, 1> Vector;
typedef Eigen::Matrix<Float,  2, 1> Vector2;
typedef Eigen::Matrix<float,  2, 1> Vector2f;
typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef Eigen::Matrix<int,    2, 1> Vector2i;

typedef Eigen::Matrix<Float,  3, 1> Vector3;
typedef Eigen::Matrix<float,  3, 1> Vector3f;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<int,    3, 1> Vector3i;

typedef Eigen::Matrix<Float,  4, 1> Vector4;
typedef Eigen::Matrix<float,  4, 1> Vector4f;
typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<int,    4, 1> Vector4i;

typedef Eigen::Matrix<Float,  2, 2> Matrix2x2;
typedef Eigen::Matrix<Float,  3, 3> Matrix3x3;
typedef Eigen::Matrix<Float,  4, 4> Matrix4x4;

typedef Eigen::Matrix<Float,  Eigen::Dynamic, Eigen::Dynamic> DynamicMatrix;
typedef Eigen::Matrix<Float,  Eigen::Dynamic, 1> DynamicVector;

template <typename T, int N> class TPoint;

typedef TPoint<Float,  3> Point;
typedef TPoint<Float,  2> Point2;
typedef TPoint<float,  2> Point2f;
typedef TPoint<double, 2> Point2d;
typedef TPoint<int,    2> Point2i;
typedef TPoint<Float,  3> Point3;
typedef TPoint<float,  3> Point3f;
typedef TPoint<double, 3> Point3d;
typedef TPoint<int,    3> Point3i;
typedef TPoint<Float,  4> Point4;
typedef TPoint<float,  4> Point4f;
typedef TPoint<double, 4> Point4d;
typedef TPoint<int,    4> Point4i;

MTS_NAMESPACE_END

#endif /* __VECTOR_H */
