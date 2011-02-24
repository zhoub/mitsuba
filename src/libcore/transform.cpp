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

#include <mitsuba/core/transform.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Linear transformation class
// -----------------------------------------------------------------------

Transform Transform::operator*(const Transform &t) const {
	return Transform(m_transform * t.m_transform,
		t.m_invTransform * m_invTransform);
}

Transform Transform::translate(const Vector &v) {
	Matrix4x4 trafo, invTrafo;
	
	trafo <<
		1, 0, 0, v.x(),
		0, 1, 0, v.y(),
		0, 0, 1, v.z(),
		0, 0, 0, 1;
	
	invTrafo <<
		1, 0, 0, -v.x(),
		0, 1, 0, -v.y(),
		0, 0, 1, -v.z(),
		0, 0, 0, 1;

	return Transform(trafo, invTrafo);
}

Transform Transform::scale(const Vector &v) {
	Matrix4x4 trafo, invTrafo;

	trafo <<
		v.x(), 0,     0,     0,
		0,     v.y(), 0,     0,
		0,     0,     v.z(), 0,
		0,     0,     0,     1;

	invTrafo <<
		1.0f/v.x(), 0,          0,          0,
		0,          1.0f/v.y(), 0,          0,
		0,          0,          1.0f/v.z(), 0,
		0,          0,          0,          1;

	return Transform(trafo, invTrafo);
}

Transform Transform::rotate(const Vector &axis, Float angle) {
	/* Make sure that the axis is normalized */
	Vector naxis = axis.normalized();

	Float rotAngle = degToRad(angle);
	Float sinAngle = std::sin(rotAngle);
	Float cosAngle = std::cos(rotAngle);

	Matrix4x4 result;

	result << naxis.x() * naxis.x() + (1.f - naxis.x() * naxis.x()) * cosAngle,
			  naxis.x() * naxis.y() * (1.f - cosAngle) - naxis.z() * sinAngle,
			  naxis.x() * naxis.z() * (1.f - cosAngle) + naxis.y() * sinAngle,
			  0,
			  naxis.x() * naxis.y() * (1.f - cosAngle) + naxis.z() * sinAngle,
			  naxis.y() * naxis.y() + (1.f - naxis.y() * naxis.y()) * cosAngle,
			  naxis.y() * naxis.z() * (1.f - cosAngle) - naxis.x() * sinAngle,
			  0,
			  naxis.x() * naxis.z() * (1.f - cosAngle) - naxis.y() * sinAngle,
			  naxis.y() * naxis.z() * (1.f - cosAngle) + naxis.x() * sinAngle,
			  naxis.z() * naxis.z() + (1.f - naxis.z() * naxis.z()) * cosAngle,
			  0,
			  0, 0, 0, 1;

	return Transform(result, result.transpose());
}

Transform Transform::perspective(Float fov, Float clipNear, Float clipFar) {
	/* Project vectors in camera space onto a plane at z=1:
	 *
	 *  xProj = x / z
	 *  yProj = y / z
	 *  zProj = (far * (z - near)) / (z * (far-near))
	 *  
	 *  Camera-space depths are not mapped linearly!
	 */

	Float recip = 1.0f / (clipFar - clipNear);

	Matrix4x4 trafo;

	trafo <<
		1,   0,   0,   0,
		0,   1,   0,   0,
		0,   0,   clipFar * recip, -clipNear * clipFar * recip,
		0,   0,   1,   0;

	/* Perform a scale so that the field of view is mapped
	 * to the interval [-1, 1] */
	Float cot = 1.0f / std::tan(degToRad(fov / 2.0f));

	return Transform::scale(Vector(cot, cot, 1.0f)) * Transform(trafo);
}

Transform Transform::glPerspective(Float fov, Float clipNear, Float clipFar) {
	Float recip = 1.0f / (clipNear - clipFar);
	Float cot = 1.0f / std::tan(degToRad(fov / 2.0f));

	Matrix4x4 trafo;

	trafo <<
		cot,   0,     0,   0,
		0,     cot,   0,   0,
		0,     0,     (clipFar + clipNear) * recip,  2 * clipFar * clipNear * recip,
		0,     0,     -1,   0;

	return Transform(trafo);
}

Transform Transform::glFrustum(Float left, Float right, Float bottom, Float top, Float nearVal, Float farVal) {
	Float invFMN = 1 / (farVal-nearVal);
	Float invTMB = 1 / (top-bottom);
	Float invRML = 1 / (right-left);

	Matrix4x4 trafo;
	
	trafo <<
		2*nearVal*invRML, 0, (right+left)*invRML, 0,
		0, 2*nearVal*invTMB, (top+bottom)*invTMB, 0,
		0, 0, -(farVal + nearVal) * invFMN, -2*farVal*nearVal*invFMN,
		0, 0, -1, 0;

	return Transform(trafo);
}

Transform Transform::orthographic(Float clipNear, Float clipFar) {
	return scale(Vector(1.0f, 1.0f, 1.0f / (clipFar - clipNear))) *
		   translate(Vector(0.0f, 0.0f, -clipNear));
}

Transform Transform::glOrthographic(Float clipNear, Float clipFar) {
	Float a = -2.0f / (clipFar - clipNear),
	      b = (clipFar + clipNear) / (clipFar - clipNear);

	Matrix4x4 trafo;
	
	trafo <<
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, a, b,
		0, 0, 0, 1;

	return Transform(trafo);
}


Transform Transform::lookAt(const Point &p, const Point &t, const Vector &up) {
	Matrix4x4 result;

	Vector dir = (t-p).normalized();
	Vector right = dir.cross(up).normalized();

	/* Generate a new, orthogonalized up vector */
	Vector newUp = right.cross(dir);

	result <<
		right, newUp, dir, p,
		0, 0, 0, 1;

	return Transform(result);
}

Transform Transform::fromFrame(const Frame &frame) {
	Matrix4x4 result;

	result <<
		frame.s, frame.t, frame.n, Vector::Zero(),
		0, 0, 0, 1;

	return Transform(result, result.transpose());
}

MTS_NAMESPACE_END
