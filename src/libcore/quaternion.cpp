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

#include <mitsuba/core/quaternion.h>

MTS_NAMESPACE_BEGIN

Quaternion Quaternion::exp() const {
	Float theta = v.norm();
	Float c = std::cos(theta);

	if (theta > Epsilon) 
		return Quaternion(v * (std::sin(theta) / theta), c);
	else
		return Quaternion(v, c);
}

Quaternion Quaternion::log() const {
	Float scale = v.norm();
	Float theta = std::atan2(scale, w);

	if (scale > 0)
		scale = theta/scale;

	return Quaternion(v * scale, 0.0f);
}

Quaternion Quaternion::fromAxisAngle(const Vector &axis, Float angle) {
	Float sinValue = std::sin(angle/2.0f), cosValue = std::cos(angle/2.0f);
	return Quaternion(axis.normalized() * sinValue, cosValue);
}

Quaternion Quaternion::fromDirectionPair(const Vector &from, const Vector &to) {
	Float dp = from.dot(to);
	if (dp > 1-Epsilon) {
		// there is nothing to do
		return Quaternion();
	} else if (dp < -(1-Epsilon)) {
		// Use a better-conditioned method for opposite directions
		Vector rotAxis = from.cross(Vector::UnitX());
		Float length = rotAxis.norm();
		if (length < Epsilon) {
			rotAxis = from.cross(Vector::UnitY());
			length = rotAxis.norm();
		}
		rotAxis /= length;
		return Quaternion(rotAxis, 0);
	} else {
		// Find cos(theta) and sin(theta) using half-angle formulae
		Float cosTheta = std::sqrt(0.5f * (1 + dp));
		Float sinTheta = std::sqrt(0.5f * (1 - dp));
		Vector rotAxis = from.cross(to).normalized();
		return Quaternion(rotAxis * sinTheta, cosTheta);
	}
}

Quaternion Quaternion::fromTransform(const Transform trafo) {
	/// Implementation from PBRT
	const Matrix4x4 &m = trafo.getMatrix();
	Float trace = m.trace();
	Vector v; Float w;
	if (trace > 0.f) {
		// Compute w from the matrix trace, then xyz
		// 4w^2 = m[0, 0] + m[1, 1] + m[2, 2] + m[3, 3] (but m[3, 3] == 1)
		Float s = std::sqrt(trace + 1.0f);
		w = s / 2.0f;
		s = 0.5f / s;
		v.x() = (m(2, 1) - m(1, 2)) * s;
		v.y() = (m(0, 2) - m(2, 0)) * s;
		v.z() = (m(1, 0) - m(0, 1)) * s;
	} else {
		// Compute largest of $x$, $y$, or $z$, then remaining components
		const int nxt[3] = {1, 2, 0};
		Float q[3];
		int i = 0;
		if (m(1, 1) > m(0, 0)) i = 1;
		if (m(2, 2) > m(i, i)) i = 2;
		int j = nxt[i];
		int k = nxt[j];
		Float s = std::sqrt((m(i, i) - (m(j, j) + m(k, k))) + 1.0);
		q[i] = s * 0.5f;
		if (s != 0.f) s = 0.5f / s;
		w = (m(k, j) - m(j, k)) * s;
		q[j] = (m(j, i) + m(i, j)) * s;
		q[k] = (m(k, i) + m(i, k)) * s;
		v.x() = q[0];
		v.y() = q[1];
		v.z() = q[2];
	}
	return Quaternion(v, w);
}

Quaternion Quaternion::fromEulerAngles(EEulerAngleConvention conv,
		Float x, Float y, Float z) {
	Quaternion qx = fromAxisAngle(Vector::UnitX(), x);
	Quaternion qy = fromAxisAngle(Vector::UnitY(), y);
	Quaternion qz = fromAxisAngle(Vector::UnitZ(), z);

	switch (conv) {
		case EEulerXYZ:
			return qz * qy * qx;
		case EEulerXZY:
			return qy * qz * qx;
		case EEulerYXZ:
			return qz * qx * qy;
		case EEulerYZX:
			return qx * qz * qy;
		case EEulerZXY:
			return qy * qx * qz;
		case EEulerZYX:
			return qx * qy * qz;
		default:
			SLog(EError, "Internal error!");
			return Quaternion();
	}
}

Transform Quaternion::toTransform() const {
	/// Implementation from PBRT
	Float xx = v.x() * v.x(), yy = v.y() * v.y(), zz = v.z() * v.z();
	Float xy = v.x() * v.y(), xz = v.x() * v.z(), yz = v.y() * v.z();
	Float wx = v.x() * w,   wy = v.y() * w,   wz = v.z() * w;

	Matrix4x4 m;

	m <<
		1.0f - 2.0f * (yy + zz),
		2.0f * (xy + wz),
		2.0f * (xz - wy),
		0.0f,

		2.0f * (xy - wz),
		1.0f - 2.f * (xx + zz),
		2.0f * (yz + wx),
		0.0f,

		2.f * (xz + wy),
		2.f * (yz - wx),
		1.f - 2.f * (xx + yy),
		0.0f,

		0.0f, 0.0f, 0.0f, 1.0f;

	return Transform(m.transpose(), m);
}
	
Quaternion Quaternion::slerp(const Quaternion &q, Float t) const {
	Float cosTheta = dot(q);
	if (cosTheta > .9995f) {
		// Revert to plain linear interpolation
		return (*this * (1.0f - t) +  q * t).normalized();
	} else {
		Float theta = std::acos(clamp(cosTheta, -1.0f, 1.0f));
		Float thetap = theta * t;
		Quaternion qperp = (q - *this * cosTheta).normalized();
		return *this * std::cos(thetap) + qperp * std::sin(thetap);
	}
}

MTS_NAMESPACE_END
