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

#if !defined(__TRANSFORM_H)
#define __TRANSFORM_H

#include <mitsuba/core/ray.h>
#include <Eigen/LU>

MTS_NAMESPACE_BEGIN

/**
 * \brief Encapsulates a 4x4 transformation in homogeneous coordinates
 * and its inverse.
 *
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE Transform {
public:
	/// Create an identity transformation
	inline Transform() : m_transform(Matrix4x4::Identity()),
			m_invTransform(Matrix4x4::Identity()) { }

	/// Unserialize a transformation from a stream
	inline Transform(Stream *stream) { 
		m_transform = Matrix4x4(stream);
		m_invTransform = Matrix4x4(stream);
	}

	/** \brief Create a transformation from the given matrix
	 * and internally compute the inverse
	 */
	explicit Transform(const Matrix4x4 &trafo)
		: m_transform(trafo) {
		bool success = false;
		trafo.computeInverseWithCheck(m_invTransform, success);
		if (!success)
			SLog(EError, "Unable to invert singular matrix %s", trafo.toString().c_str());
	}

	/// Create a transformation from the given matrices
	inline Transform(const Matrix4x4 &trafo, const Matrix4x4 &invTrafo)
		: m_transform(trafo), m_invTransform(invTrafo) { }

	/// Return the inverse transformation
	Transform inverse() const {
		return Transform(m_invTransform, m_transform);
	}

	/// Matrix-matrix multiplication
	Transform operator*(const Transform &t) const;

	/// Return the determinant of the upper left 3x3 submatrix
	inline Float det3x3() const {
		return m_transform.topLeftCorner<3, 3>().determinant();
	}

	/// Test for a scale component
	inline bool hasScale() const {
		for (int i=0; i<3; ++i) {
			for (int j=i; j<3; ++j) {
				Float dp = 
					m_transform.row(i).head<3>().dot(
					m_transform.col(j ).head<3>());

				if (i == j && std::abs(dp-1) > Epsilon)
					return true;
				else if (i != j && std::abs(dp) > Epsilon)
					return true;
			}
		}
		return false;
	}

	/// Test if this is the identity matrix
	inline bool isIdentity() const {
		return m_transform == Matrix4x4::Identity();
	}

	/// Transform a point by an arbitrary matrix in homogeneous coordinates
	inline Point operator()(const Point &p) const {
		Vector4 result = m_transform * Vector4(p[0], p[1], p[2], (Float) 1.0f);

#ifdef MTS_DEBUG
		if (result.w() == 0)
			SLog(EWarn, "w==0 in Transform::operator(Point &)");
#endif

		if (result.w() == 1.0f)
			return result.head<3>();
		else
			return result.head<3>() / result.w();
	}

	/// Transform a point by an affine / non-projective matrix
	inline Point transformAffine(const Point &p) const {
		return (m_transform * Vector4(p[0], p[1], p[2], (Float) 1.0f)).head<3>();
	}

	/// Transform a point by an arbitrary matrix in homogeneous coordinates
    inline void operator()(const Point &p, Point &dest) const {
		Vector4 result = m_transform * Vector4(p[0], p[1], p[2], (Float) 1.0f);

#ifdef MTS_DEBUG
		if (result.w() == 0)
			SLog(EWarn, "w==0 in Transform::operator(Point &, Point &)");
#endif
		if (result.w() != 1.0f)
			dest = result.head<3>() / result.w();
		else
			dest = result.head<3>();
	}

	/// Transform a vector by the upper left 3x3 submatrix
    inline Vector operator()(const Vector &v) const {
		return m_transform.topLeftCorner<3, 3>() * v;
	}

	/// Transform a vector by the upper left 3x3 submatrix
    inline void operator()(const Vector &v, Vector &dest) const {
		dest = m_transform.topLeftCorner<3, 3>() * v;
	}

	/// Transform a normal vector by the inverse transpose matrix
    inline Normal operator()(const Normal &n) const {
		return m_invTransform.topLeftCorner<3, 3>().transpose() * n;
	}

	/// Transform a normal vector by the inverse transpose matrix
    inline void operator()(const Normal &n, Normal &dest) const {
		dest = m_invTransform.topLeftCorner<3, 3>().transpose() * n;
	}

	/// 4D matrix-vector multiplication
	inline Vector4 operator()(const Vector4 &v) const {
		return m_transform * v;
	}

	/// 4D matrix-vector multiplication
	inline void operator()(const Vector4 &v, Vector4 &dest) const {
		dest = m_transform * v;
	}

	/// Transform a ray
	inline void operator()(const Ray &a, Ray &b) const {
		b.mint = a.mint;
		b.maxt = a.maxt;
		operator()(a.o, b.o);
		operator()(a.d, b.d);
#ifdef MTS_DEBUG_FP
		bool state = disableFPExceptions();
#endif
		/* Re-compute the reciprocal */
		b.dRcp = b.d.cwiseInverse();
		b.time = a.time;
#ifdef MTS_DEBUG_FP
		restoreFPExceptions(state);
#endif
	}
	
	/// Return the underlying matrix
	inline const Matrix4x4 &getMatrix() const { return m_transform; }

	/// Return the underlying inverse matrix (const version)
	inline const Matrix4x4 &getInverseMatrix() const { return m_invTransform; }

	/// Create a translation transformation
	static Transform translate(const Vector &v);

	/// Create a rotation transformation around an arbitrary axis. The angle is specified in degrees
	static Transform rotate(const Vector &axis, Float angle);

	/// Create a scale transformation
	static Transform scale(const Vector &v);
	
	/** \brief Create a perspective transformation.
	 *   (Maps [near, far] to [0, 1])
	 * \param fov Field of view in degrees
	 * \param clipNear Near clipping plane
	 * \param clipFar Far clipping plane
	 */
	static Transform perspective(Float fov, Float clipNear, Float clipFar);
	
	/** \brief Create a perspective transformation for OpenGL.
	 *   (Maps [near, far] to [-1, 1])
	 * \param fov Field of view in degrees
	 * \param clipNear Near clipping plane distance
	 * \param clipFar Far clipping plane distance
	 */
	static Transform glPerspective(Float fov, Float clipNear, Float clipFar);

	/** \brief Create a perspective transformation for OpenGL.
	 * \param left Left clipping plane coordinate
	 * \param right Right clipping plane coordinate
	 * \param top Top clipping plane coordinate
	 * \param bottom Bottom clipping plane coordinate
	 * \param nearVal Near clipping plane distance
	 * \param farVal Far clipping plane distance
	 */
	static Transform glFrustum(Float left, Float right, Float bottom, Float top, Float nearVal, Float farVal);

	/** \brief Create an orthographic transformation, which maps Z to [0,1]
	 * and leaves the X and Y coordinates untouched.
	 * \param clipNear Near clipping plane
	 * \param clipFar Far clipping plane
	 */
	static Transform orthographic(Float clipNear, Float clipFar);
	
	/** \brief Create an orthographic transformation for OpenGL
	 * \param clipNear Near clipping plane
	 * \param clipFar Far clipping plane
	 */
	static Transform glOrthographic(Float clipNear, Float clipFar);

	/** \brief Create a look-at camera transformation
	 * \param p Camera position
	 * \param t Target vector
	 * \param u Up vector
	 */
	static Transform lookAt(const Point &p, const Point &t, const Vector &u);

	/** \brief Create an orthogonal transformation that takes
	 * the standard to the supplied frame
	 */
	static Transform fromFrame(const Frame &frame);

	/// Serialize a transformation to a stream
	inline void serialize(Stream *stream) const {
		m_transform.serialize(stream);
		m_invTransform.serialize(stream);
	}

	/// Return a string representation
	std::string toString() const {
		return m_transform.toString();
	}
private:
	Matrix4x4 m_transform;
	Matrix4x4 m_invTransform;
};

MTS_NAMESPACE_END

#endif /* __TRANSFORM_H */
