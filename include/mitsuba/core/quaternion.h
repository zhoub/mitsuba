/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUFloat ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__QUATERNION_H)
#define __QUATERNION_H

#include <mitsuba/core/transform.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Quaternion data structure
 * \ingroup libcore
 */
struct Quaternion {
	/// Used by \ref Quaternion::fromEulerAngles
	enum EEulerAngleConvention {
		EEulerXYZ = 0,
		EEulerXZY,
		EEulerYXZ,
		EEulerYZX,
		EEulerZXY,
		EEulerZYX
	};

	/// Imaginary component
	Vector v;

	/// Real component
	Float w;

	/// Create a unit quaternion
	inline Quaternion() : v(Vector::Zero()), w(1) { }

	/**
	 * Initialize the quaternion with the specified 
	 * real and imaginary components
	 */
	inline Quaternion(const Vector &v, Float w) : v(v), w(w) {  }

	/// Unserialize a quaternion from a binary data stream
	explicit Quaternion(Stream *stream) {
		v = Vector(stream);
		w = stream->readFloat();
	}

	/// Add two quaternions and return the result
	inline Quaternion operator+(const Quaternion &q) const {
		return Quaternion(v + q.v, w + q.w);
	}

	/// Subtract two quaternions and return the result
	inline Quaternion operator-(const Quaternion &q) const {
		return Quaternion(v - q.v, w - q.w);
	}

	/// Add a quaternion to the current one
	inline Quaternion& operator+=(const Quaternion &q) {
		v += q.v; w += q.w; 
		return *this;
	}

	/// Subtract a quaternion
	inline Quaternion& operator-=(const Quaternion &q) {
		v -= q.v; w -= q.w;
		return *this;
	}

	/// Multiply the quaternion by the given scalar and return the result
	inline Quaternion operator*(Float f) const {
		return Quaternion(v*f, w*f);
	}

	/// Multiply the quaternion by the given scalar
	inline Quaternion &operator*=(Float f) {
		v *= f; w *= f; 
		return *this;
	}

	/// Divide the quaternion by the given scalar and return the result
	inline Quaternion operator/(Float f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Quaternion: Division by zero!");
#endif
		Float recip = (Float) 1 / f;
		return Quaternion(v * recip, w * recip);
	}

	/// Divide the quaternion by the given scalar
	inline Quaternion &operator/=(Float f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Quaternion: Division by zero!");
#endif
		Float recip = (Float) 1 / f;
		v *= recip; w *= recip; 
		return *this;
	}

	/// Quaternion multiplication
	inline Quaternion &operator*=(const Quaternion &q) {
		Float tmp = w * q.w - v.dot(q.v);
		v = v.cross(q.v) + q.w * v + w * q.v;
		w = tmp;
		return *this;
	}

	/// Quaternion multiplication (creates a temporary)
	Quaternion operator*(const Quaternion &q) const {
		return Quaternion(v.cross(q.v) + q.w * v + w * q.v,
			w * q.w - v.dot(q.v));
	}

	/// Equality test
	bool operator==(const Quaternion &q) const {
		return v == q.v && w == q.w;
	}

	/// Inequality test
	bool operator!=(const Quaternion &q) const {
		return v != q.v || w != q.w;
	}

	inline Float dot(const Quaternion &q) const {
		return v.dot(q.v) + w * q.w;
	}

	inline void normalize() {
		*this /= dot(*this);
	}

	inline Quaternion normalized() const {
		return *this / dot(*this);
	}


	/// Return the rotation axis of this quaternion
    inline Vector getAxis() const { return v.normalized(); }

	/// Return the rotation angle of this quaternion (in radians)
    inline Float getAngle() const { return 2 * std::acos(w); }

	/**
	 * \brief Compute the exponential of a quaternion with
	 * scalar part w = 0. 
	 *
	 * Based on code the appendix of
	 * "Quaternion Calculus for Computer Graphics" by Ken Shoemake
	 */
	Quaternion exp() const;

	/**
	 * \brief Compute the natural logarithm of a unit quaternion
	 *
	 * Based on code the appendix of
	 * "Quaternion Calculus for Computer Graphics" by Ken Shoemake
	 */
	Quaternion log() const;
	
	/**
	 * \brief Construct an unit quaternion, which represents a rotation
	 * around \a axis by \a angle radians.
	 */
	static Quaternion fromAxisAngle(const Vector &axis, Float angle);

	/**
	 * \brief Construct an unit quaternion, which rotates unit direction 
	 * \a from onto \a to.
	 */
	static Quaternion fromDirectionPair(const Vector &from, const Vector &to);
	
	/**
	 * \brief Construct an unit quaternion matching the supplied
	 * rotation matrix.
	 */
	static Quaternion fromTransform(const Transform trafo);
	/**
	 * \brief Construct an unit quaternion matching the supplied
	 * rotation expressed in Euler angles (in radians)
	 */
	static Quaternion fromEulerAngles(EEulerAngleConvention conv,
			Float x, Float y, Float z);

	/// Compute the rotation matrix for the given quaternion
	Transform toTransform() const;

	/// Spherlical linear interpolation, where q is in [0, 1]
	Quaternion slerp(const Quaternion &q, Float t) const;

	/// Serialize this quaternion to a binary data stream
	void serialize(Stream *stream) const {
		v.serialize(stream);
		stream->writeFloat(w);
	}

	/// Return a readable string representation of this quaternion
	std::string toString() const {
		std::ostringstream oss;
		oss << "Quaternion[v=" << v.toString() << ", w=" << w << "]";
		return oss.str();
	}
};

inline Quaternion operator*(Float f, const Quaternion &v) {
	return v*f;
}

MTS_NAMESPACE_END

#endif /* __QUATERNION_H */
