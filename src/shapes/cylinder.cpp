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

#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

class Cylinder : public Shape {
private:
	Transform m_objectToWorld;
	Transform m_worldToObject;
	Float m_radius, m_length, m_invSurfaceArea;
public:
	Cylinder(const Properties &props) : Shape(props) {
		/**
		 * There are two ways of instantiating cylinders: either,
		 * one can specify a linear transformation to from the
		 * unit cylinder using the 'toWorld' parameter, or one
		 * can explicitly specify two points and a radius.
		 */
		if (props.hasProperty("p1") && props.hasProperty("p2")
				&& props.hasProperty("radius")) {
			Point p1 = props.getPoint("p1"), p2 = props.getPoint("p2");
			Vector rel = p2 - p1;
			Float radius = props.getFloat("radius");
			Float length = rel.norm();

			m_objectToWorld = 
				Transform::translate(Vector(p1)) *
				Transform::fromFrame(Frame(rel/length));
			m_radius = radius;
			m_length = length;
		} else {
			Transform objectToWorld = props.getTransform("toWorld", Transform());
			m_radius = objectToWorld(Vector(1,0,0)).norm();
			m_length = objectToWorld(Vector(0,0,1)).norm();
			// Remove the scale from the object-to-world trasnsform
			m_objectToWorld = objectToWorld * Transform::scale(
					Vector(1/m_radius, 1/m_radius, 1/m_length));
		}
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
		Assert(m_length > 0 && m_radius > 0);
	}

	Cylinder(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		m_objectToWorld = Transform(stream);
		m_radius = stream->readFloat();
		m_length = stream->readFloat();
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
		stream->writeFloat(m_radius);
		stream->writeFloat(m_length);
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		Ray ray;

		/* Transform into the local coordinate system and normalize */
		m_worldToObject(_ray, ray);

		const Float
			ox = ray.o.x,
			oy = ray.o.y,
			dx = ray.d.x, 
			dy = ray.d.y;

		const Float A = dx*dx + dy*dy;
		const Float B = 2 * (dx*ox + dy*oy);
		const Float C = ox*ox + oy*oy - m_radius*m_radius;

		Float nearT, farT;
		if (!solveQuadratic(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;

		const Float zPosNear = ray.o.z + ray.d.z * nearT;
		const Float zPosFar = ray.o.z + ray.d.z * farT;
		if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
			t = nearT;
		} else if (zPosFar >= 0 && zPosFar <= m_length) {
			if (farT > maxt)
				return false;
			t = farT;
		} else {
			return false;
		}

		return true;
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
		Ray ray;

		/* Transform into the local coordinate system and normalize */
		m_worldToObject(_ray, ray);

		const Float
			ox = ray.o.x,
			oy = ray.o.y,
			dx = ray.d.x, 
			dy = ray.d.y;

		const Float A = dx*dx + dy*dy;
		const Float B = 2 * (dx*ox + dy*oy);
		const Float C = ox*ox + oy*oy - m_radius*m_radius;

		Float nearT, farT;
		if (!solveQuadratic(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;

		const Float zPosNear = ray.o.z + ray.d.z * nearT;
		const Float zPosFar = ray.o.z + ray.d.z * farT;
		if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
			return true;
		} else if (zPosFar >= 0 && zPosFar <= m_length && farT <= maxt) {
			return true;
		} else {
			return false;
		}
	}

	void fillIntersectionRecord(const Ray &ray,
			const void *temp, Intersection &its) const {
		its.p = ray(its.t);

		Point local = m_worldToObject(its.p);
		Float phi = std::atan2(local.y, local.x);
		if (phi < 0)
			phi += 2*M_PI;
		its.uv.x = local.z / m_length;
		its.uv.y = phi / (2*M_PI);

		Vector dpdu = Vector(-local.y, local.x, 0) * (2*M_PI);
		Vector dpdv = Vector(0, 0, m_length);
		its.dpdu = m_objectToWorld(dpdu);
		its.dpdv = m_objectToWorld(dpdv);
		its.geoFrame.n = Normal(m_objectToWorld(dpdu.cross(dpdv)).normalized());
		its.geoFrame.s = its.dpdu.normalized();
		its.geoFrame.t = its.dpdv.normalized();
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.hasUVPartials = false;
		its.shape = this;
	}

	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
		Point p = Point(m_radius * std::cos(sample.y), 
			m_radius * std::sin(sample.y), 
			sample.x * m_length);
		sRec.p = m_objectToWorld(p);
		sRec.n = m_objectToWorld(Normal(p.x, p.y, 0.0f)).normalized();
		return m_invSurfaceArea;
	}

	inline BoundingBox3 getBoundingBox() const {
		Vector x1 = m_objectToWorld(Vector(m_radius, 0, 0));
		Vector x2 = m_objectToWorld(Vector(0, m_radius, 0));
		Point p0 = m_objectToWorld(Point::Zero());
		Point p1 = m_objectToWorld(Point(0, 0, m_length));
		BoundingBox3 result;

		/* To bound the cylinder, it is sufficient to find the
		   smallest box containing the two circles at the endpoints.
		   This can be done component-wise as follows */

		for (int i=0; i<3; ++i) {
			Float range = std::sqrt(x1[i]*x1[i] + x2[i]*x2[i]);

			result.min[i] = std::min(std::min(result.min[i], 
						p0[i]-range), p1[i]-range);
			result.max[i] = std::max(std::max(result.max[i], 
						p0[i]+range), p1[i]+range);
		}

		return result;
	}

	/**
	 * Compute the ellipse created by the intersection of an infinite
	 * cylinder and a plane. Returns false in the degenerate case.
	 * Based on:
	 * www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
	 */
	bool intersectCylPlane(Point planePt, Normal planeNrml,
			Point cylPt, Vector cylD, Float radius, Point &center,
			Vector *axes, Float *lengths) const {
		if (std::abs(planeNrml.dot(cylD)) < Epsilon)
			return false;

		Vector B, A = cylD - cylD.dot(planeNrml)*planeNrml;

		Float length = A.norm();
		if (length != 0) {
			A /= length;
			B = cross(planeNrml, A);
		} else {
			coordinateSystem(planeNrml, A, B);
		}

		Vector delta = planePt - cylPt,
			   deltaProj = delta - cylD*delta.dot(cylD);

		Float aDotD = A.dot(cylD);
		Float bDotD = B.dot(cylD);
		Float c0 = 1-aDotD*aDotD;
		Float c1 = 1-bDotD*bDotD;
		Float c2 = 2*A.dot(deltaProj);
		Float c3 = 2*B.dot(deltaProj);
		Float c4 = delta.dot(deltaProj) - radius*radius;

		Float lambda = (c2*c2/(4*c0) + c3*c3/(4*c1) - c4)/(c0*c1);

		Float alpha0 = -c2/(2*c0),
			  beta0 = -c3/(2*c1);

		lengths[0] = std::sqrt(c1*lambda),
		lengths[1] = std::sqrt(c0*lambda);

		center = planePt + alpha0 * A + beta0 * B;
		axes[0] = A;
		axes[1] = B;
		return true;
	}

	BoundingBox3 intersectCylFace(int axis,
			const Point &min, const Point &max,
			const Point &cylPt, const Vector &cylD) const {
		int axis1 = (axis + 1) % 3;
		int axis2 = (axis + 2) % 3;

		Normal planeNrml(Normal::Zero());
		planeNrml[axis] = 1;

		Point ellipseCenter;
		Vector ellipseAxes[2];
		Float ellipseLengths[2];

		BoundingBox3 bbox;
		if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius, 
			ellipseCenter, ellipseAxes, ellipseLengths)) {
			/* Degenerate case -- return an invalid bounding box. This is
			   not a problem, since one of the other faces will provide
			   enough information to arrive at a correct clipped bounding box */
			return bbox;
		}

		/* Intersect the ellipse against the sides of the bounding box face */
		for (int i=0; i<4; ++i) {
			Point p1, p2;
			p1[axis] = p2[axis] = min[axis];
			p1[axis1] = ((i+1) & 2) ? min[axis1] : max[axis1];
			p1[axis2] = ((i+0) & 2) ? min[axis2] : max[axis2];
			p2[axis1] = ((i+2) & 2) ? min[axis1] : max[axis1];
			p2[axis2] = ((i+1) & 2) ? min[axis2] : max[axis2];

			Point2 p1l(
				(p1 - ellipseCenter).dot(ellipseAxes[0]) / ellipseLengths[0],
				(p1 - ellipseCenter).dot(ellipseAxes[1]) / ellipseLengths[1]);
			Point2 p2l(
				(p2 - ellipseCenter).dot(ellipseAxes[0]) / ellipseLengths[0],
				(p2 - ellipseCenter).dot(ellipseAxes[1]) / ellipseLengths[1]);

			Vector2 rel = p2l-p1l;
			Float A = rel.squaredNorm();
			Float B = 2*p1l.dot(rel);
			Float C = p1l.squaredNorm()-1;

			Float x0, x1;
			if (solveQuadratic(A, B, C, x0, x1)) {
				if (x0 >= 0 && x0 <= 1)
					bbox.expandBy(p1+(p2-p1)*x0);
				if (x1 >= 0 && x1 <= 1)
					bbox.expandBy(p1+(p2-p1)*x1);
			}
		}

		ellipseAxes[0] *= ellipseLengths[0];
		ellipseAxes[1] *= ellipseLengths[1];
		BoundingBox3 faceBounds(min, max);

		/* Find the componentwise maxima of the ellipse */
		for (int i=0; i<2; ++i) {
			int j = (i==0) ? axis1 : axis2;
			Float alpha = ellipseAxes[0][j];
			Float beta = ellipseAxes[1][j];
			Float ratio = beta/alpha, tmp = std::sqrt(1+ratio*ratio);
			Float cosTheta = 1/tmp, sinTheta = ratio/tmp;
			Point p1 = ellipseCenter + cosTheta*ellipseAxes[0] + sinTheta*ellipseAxes[1];
			Point p2 = ellipseCenter - cosTheta*ellipseAxes[0] - sinTheta*ellipseAxes[1];

			if (faceBounds.contains(p1)) 
				bbox.expandBy(p1);
			if (faceBounds.contains(p2)) 
				bbox.expandBy(p2);
		}

		return bbox;
	}

	BoundingBox3 getClippedBoundingBox(const BoundingBox3 &box) const {
		/* Compute a base bounding box */
		BoundingBox3 base(getBoundingBox());
		base.clip(box);
		
		Point cylPt = m_objectToWorld(Point::Zero());
		Vector cylD(m_objectToWorld(Vector::UnitZ()));

		/* Now forget about the cylinder ends and 
		   intersect an infinite cylinder with each bounding box face */
		BoundingBox3 clippedBoundingBox;
		clippedBoundingBox.expandBy(intersectCylFace(0, 
				Point(base.min.x(), base.min.y(), base.min.z()),
				Point(base.min.x(), base.max.y(), base.max.z()),
				cylPt, cylD));

		clippedBoundingBox.expandBy(intersectCylFace(0,
				Point(base.max.x(), base.min.y(), base.min.z()),
				Point(base.max.x(), base.max.y(), base.max.z()),
				cylPt, cylD));

		clippedBoundingBox.expandBy(intersectCylFace(1, 
				Point(base.min.x(), base.min.y(), base.min.z()),
				Point(base.max.x(), base.min.y(), base.max.z()),
				cylPt, cylD));

		clippedBoundingBox.expandBy(intersectCylFace(1,
				Point(base.min.x(), base.max.y(), base.min.z()),
				Point(base.max.x(), base.max.y(), base.max.z()),
				cylPt, cylD));

		clippedBoundingBox.expandBy(intersectCylFace(2, 
				Point(base.min.x(), base.min.y(), base.min.z()),
				Point(base.max.x(), base.max.y(), base.min.z()),
				cylPt, cylD));

		clippedBoundingBox.expandBy(intersectCylFace(2,
				Point(base.min.x(), base.min.y(), base.max.z()),
				Point(base.max.x(), base.max.y(), base.max.z()),
				cylPt, cylD));

		clippedBoundingBox.clip(box);
		return clippedBoundingBox;
	}

	ref<TriMesh> createTriMesh() {
		/// Choice of discretization
		const size_t phiSteps = 20;
		const Float dPhi   = (2*M_PI) / phiSteps;

		ref<TriMesh> mesh = new TriMesh("Cylinder approximation",
			phiSteps*2, phiSteps*2, true, false, false);

		Point *vertices = mesh->getVertexPositions();
		Normal *normals = mesh->getVertexNormals();
		Triangle *triangles = mesh->getTriangles();
		size_t triangleIdx = 0, vertexIdx = 0;

		for (size_t phi=0; phi<phiSteps; ++phi) {
			Float sinPhi = std::sin(phi * dPhi);
			Float cosPhi = std::cos(phi * dPhi);
			int idx0 = vertexIdx, idx1 = idx0+1;
			int idx2 = (vertexIdx+2) % (2*phiSteps), idx3 = idx2+1;
			normals[vertexIdx] = m_objectToWorld(Normal(cosPhi, sinPhi, 0));
			vertices[vertexIdx++] = m_objectToWorld(Point(cosPhi*m_radius, sinPhi*m_radius, 0));
			normals[vertexIdx] = m_objectToWorld(Normal(cosPhi, sinPhi, 0));
			vertices[vertexIdx++] = m_objectToWorld(Point(cosPhi*m_radius, sinPhi*m_radius, m_length));

			triangles[triangleIdx].idx[0] = idx0;
			triangles[triangleIdx].idx[1] = idx2;
			triangles[triangleIdx].idx[2] = idx1;
			triangleIdx++;
			triangles[triangleIdx].idx[0] = idx1;
			triangles[triangleIdx].idx[1] = idx2;
			triangles[triangleIdx].idx[2] = idx3;
			triangleIdx++;
		}

		mesh->setBSDF(m_bsdf);
		mesh->setLuminaire(m_luminaire);
		mesh->configure();

		return mesh.get();
	}

#if 0
	BoundingBox3 getBoundingBox() const {
		const Point a = m_objectToWorld(Point::Zero());
		const Point b = m_objectToWorld(Point(0, 0, m_length));

		const Float r = m_radius;
		BoundingBox3 result;
		result.expandBy(a - Vector(r, r, r));
		result.expandBy(a + Vector(r, r, r));
		result.expandBy(b - Vector(r, r, r));
		result.expandBy(b + Vector(r, r, r));
		return result;
	}
#endif

	Float getSurfaceArea() const {
		return 2*M_PI*m_radius*m_length;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Cylinder[" << endl
			<< "  radius = " << m_radius << ", " << endl
			<< "  length = " << m_length << ", " << endl
			<< "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
			<< "  luminaire = " << indent(m_luminaire.toString()) << "," << endl
			<< "  subsurface = " << indent(m_subsurface.toString())
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(Cylinder, false, Shape)
MTS_EXPORT_PLUGIN(Cylinder, "Cylinder intersection primitive");
MTS_NAMESPACE_END
