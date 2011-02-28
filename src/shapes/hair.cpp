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
#include <mitsuba/render/sahkdtree3.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

#define MTS_HAIR_USE_FANCY_CLIPPING 1

MTS_NAMESPACE_BEGIN

/**
 * \brief Space-efficient acceleration structure for cylindrical hair
 * segments with miter joints. This class expects an ASCII file containing
 * a list of hairs made from segments. Each line should contain an X,
 * Y and Z coordinate separated by a space. An empty line indicates
 * the start of a new hair.
 */
class HairKDTree : public SAHKDTree3D<HairKDTree> {
	friend class GenericKDTree<BoundingBox3, SurfaceAreaHeuristic, HairKDTree>;
	friend class SAHKDTree3D<HairKDTree>;
public:
	HairKDTree(std::vector<Point> &vertices, 
			std::vector<bool> &vertexStartsFiber, Float radius)
			: m_radius(radius) {
		/* Take the supplied vertex & start fiber arrays (without copying) */
		m_vertices.swap(vertices);
		m_vertexStartsFiber.swap(vertexStartsFiber);
		m_hairCount = 0;

		/* Compute the index of the first vertex in each segment. */
		m_segIndex.reserve(m_vertices.size());
		for (size_t i=0; i<m_vertices.size()-1; i++) {
			if (m_vertexStartsFiber[i])
				m_hairCount++;
			if (!m_vertexStartsFiber[i+1])
				m_segIndex.push_back(i);
		}
		m_segmentCount = m_segIndex.size();

		Log(EDebug, "Building a kd-tree for " SIZE_T_FMT " hair vertices, "
			SIZE_T_FMT " segments, " SIZE_T_FMT " hairs", 
			m_vertices.size(), m_segmentCount, m_hairCount);

		/* Ray-cylinder intersections are expensive. Use only the
		   SAH cost as the tree subdivision stopping criterion, 
		   not the number of primitives */
		setStopPrims(0);
		setTraversalCost(10);
		setQueryCost(30);
		buildInternal();

		Log(EDebug, "Total amount of storage (kd-tree & vertex data): %s",
			memString(m_nodeCount * sizeof(KDNode) 
			+ m_indexCount * sizeof(IndexType)
			+ vertices.size() * sizeof(Point)
			+ vertexStartsFiber.size() / 8).c_str());

		/* Optimization: replace all primitive indices by the
		   associated vertex indices (this avoids an extra 
		   indirection during traversal later on) */
		for (SizeType i=0; i<m_indexCount; ++i)
			m_indices[i] = m_segIndex[m_indices[i]];

		/* Free the segIndex array, it is not needed anymore */
		std::vector<IndexType>().swap(m_segIndex);
	}

	/// Return the bounding box of the hair kd-tree
	inline const BoundingBox3 &getBoundingBox() const {
		return m_bbox;
	}

	/// Return the list of vertices underlying the hair kd-tree
	inline const std::vector<Point> &getVertices() const {
		return m_vertices;
	}

	/**
	 * Return a boolean list specifying whether a vertex 
	 * marks the beginning of a new fiber
	 */
	inline const std::vector<bool> &getStartFiber() const {
		return m_vertexStartsFiber;
	}

	/// Return the radius of the hairs stored in the kd-tree
	inline Float getRadius() const {
		return m_radius;
	}

	/// Return the total number of segments
	inline size_t getSegmentCount() const {
		return m_segmentCount;
	}
	
	/// Return the total number of hairs
	inline size_t getHairCount() const {
		return m_hairCount;
	}
	
	/// Return the total number of vertices
	inline size_t getVertexCount() const {
		return m_vertices.size();
	}

	/// Intersect a ray with all segments stored in the kd-tree
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt, 
			Float &t, void *temp) const {
		Float tempT = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;

		if (m_bbox.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<false>(ray, mint, maxt, tempT, temp)) {
					t = tempT;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * \brief Intersect a ray with all segments stored in the kd-tree
	 * (Visiblity query version)
	 */
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt) const {
		Float tempT = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;

		if (m_bbox.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL)) 
					return true;
			}
		}
		return false;
	}

#if defined(MTS_HAIR_USE_FANCY_CLIPPING)
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

		Assert(std::abs(planeNrml.norm()-1) <Epsilon);
		Vector B, A = cylD - cylD.dot(planeNrml)*planeNrml;

		Float length = A.norm();
		if (length > Epsilon && planeNrml != cylD) {
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

	/**
	 * \brief Intersect an infinite cylinder with an 
	 * bounding box face and bound the resulting clipped ellipse
	 */
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

	BoundingBox3 getBoundingBox(IndexType index) const {
		IndexType iv = m_segIndex.at(index);
		Point center;
		Vector axes[2];
		Float lengths[2];

		bool success = intersectCylPlane(firstVertex(iv), firstMiterNormal(iv), 
			firstVertex(iv), tangent(iv), m_radius, center, axes, lengths);
		Assert(success);

		BoundingBox3 result;
		axes[0] *= lengths[0]; axes[1] *= lengths[1];
		for (int i=0; i<3; ++i) {
			Float range = std::sqrt(axes[0][i]*axes[0][i] + axes[1][i]*axes[1][i]);
			result.min[i] = std::min(result.min[i], center[i]-range);
			result.max[i] = std::max(result.max[i], center[i]+range);
		}

		success = intersectCylPlane(secondVertex(iv), secondMiterNormal(iv), 
			secondVertex(iv), tangent(iv), m_radius, center, axes, lengths);
		Assert(success);

		axes[0] *= lengths[0]; axes[1] *= lengths[1];
		for (int i=0; i<3; ++i) {
			Float range = std::sqrt(axes[0][i]*axes[0][i] + axes[1][i]*axes[1][i]);
			result.min[i] = std::min(result.min[i], center[i]-range);
			result.max[i] = std::max(result.max[i], center[i]+range);
		}
		return result;
	}

	BoundingBox3 getClippedBoundingBox(IndexType index, const BoundingBox3 &box) const {
		/* Compute a base bounding box */
		BoundingBox3 base(getBoundingBox(index));
		base.clip(box);

		IndexType iv = m_segIndex.at(index);

		Point cylPt = firstVertex(iv);
		Vector cylD = tangent(iv);

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

		clippedBoundingBox.clip(base);
		return clippedBoundingBox;
	}
#else
	/// Compute the bounding box of a segment (only used during tree construction)
	BoundingBox3 getBoundingBox(int index) const {
		IndexType iv = m_segIndex.at(index);

		// cosine of steepest miter angle
		const Float cos0 = firstMiterNormal(iv).dot(tangent(iv));
		const Float cos1 = secondMiterNormal(iv).dot(tangent(iv));
		const Float maxInvCos = 1.0 / std::min(cos0, cos1);
		const Vector expandVec(m_radius * maxInvCos);

		const Point a = firstVertex(iv);
		const Point b = secondVertex(iv);

		BoundingBox3 bbox;
		bbox.expandBy(a - expandVec);
		bbox.expandBy(a + expandVec);
		bbox.expandBy(b - expandVec);
		bbox.expandBy(b + expandVec);
		return bbox;
	}

	/// Compute the clipped bounding box of a segment (only used during tree construction)
	BoundingBox3 getClippedBoundingBox(int index, const BoundingBox3 &box) const {
		BoundingBox3 bbox(getBoundingBox(index));
		bbox.clip(box);
		return bbox;
	}
#endif

	/// Return the total number of segments
	inline int getPrimitiveCount() const {
		return m_segIndex.size();
	}

	inline bool intersect(const Ray &ray, IndexType iv, 
		Float mint, Float maxt, Float &t, void *tmp) const {
		/* First compute the intersection with the infinite cylinder */
		Float nearT, farT;
		Vector axis = tangent(iv);

		// Projection of ray onto subspace normal to axis
		Vector relOrigin = ray.o - firstVertex(iv);
		Vector projOrigin = relOrigin - axis.dot(relOrigin) * axis;
		Vector projDirection = ray.d - axis.dot(ray.d) * axis;

		// Quadratic to intersect circle in projection
		const Float A = projDirection.squaredNorm();
		const Float B = 2 * projOrigin.dot(projDirection);
		const Float C = projOrigin.squaredNorm() - m_radius*m_radius;

		if (!solveQuadratic(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;

		/* Next check the intersection points against the miter planes */
		Point pointNear = ray(nearT);
		Point pointFar = ray(farT);
		if ((pointNear - firstVertex(iv)).dot(firstMiterNormal(iv)) >= 0 &&
			(pointNear - secondVertex(iv)).dot(secondMiterNormal(iv)) <= 0 &&
			nearT >= mint) {
			t = nearT;
		} else if ((pointFar - firstVertex(iv)).dot(firstMiterNormal(iv)) >= 0 &&
				(pointFar - secondVertex(iv)).dot(secondMiterNormal(iv)) <= 0) {
			if (farT > maxt)
				return false;
			t = farT;
		} else {
			return false;
		}

		IndexType *storage = static_cast<IndexType *>(tmp);
		if (storage)
			*storage = iv;

		return true;
	}
	
	inline bool intersect(const Ray &ray, IndexType iv, 
		Float mint, Float maxt) const {
		Float tempT;
		return intersect(ray, iv, mint, maxt, tempT, NULL);
	}

	/* Some utility functions */
	inline Point firstVertex(IndexType iv) const { return m_vertices[iv]; }
	inline Point secondVertex(IndexType iv) const { return m_vertices[iv+1]; }
	inline Point prevVertex(IndexType iv) const { return m_vertices[iv-1]; }
	inline Point nextVertex(IndexType iv) const { return m_vertices[iv+2]; }

	inline bool prevSegmentExists(IndexType iv) const { return !m_vertexStartsFiber[iv]; }
	inline bool nextSegmentExists(IndexType iv) const { return !m_vertexStartsFiber[iv+2]; }

	inline Vector tangent(IndexType iv) const { return (secondVertex(iv) - firstVertex(iv)).normalized(); }
	inline Vector prevTangent(IndexType iv) const { return (firstVertex(iv) - prevVertex(iv)).normalized(); }
	inline Vector nextTangent(IndexType iv) const { return (nextVertex(iv) - secondVertex(iv)).normalized(); }

	inline Vector firstMiterNormal(IndexType iv) const {
		if (prevSegmentExists(iv))
			return (prevTangent(iv) + tangent(iv)).normalized();
		else
			return tangent(iv);
	}

	inline Vector secondMiterNormal(IndexType iv) const {
		if (nextSegmentExists(iv))
			return (tangent(iv) + nextTangent(iv)).normalized();
		else
			return tangent(iv);
	}

	MTS_DECLARE_CLASS()
protected:
	std::vector<Point> m_vertices;
	std::vector<bool> m_vertexStartsFiber;
	std::vector<IndexType> m_segIndex;
	size_t m_segmentCount;
	size_t m_hairCount;
	Float m_radius;
};

class Hair : public Shape {
public:
	Hair(const Properties &props) : Shape(props) {
		fs::path path = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		Float radius = props.getFloat("radius", 0.05f);
		/* Skip segments, whose tangent differs by less than one degree
		   compared to the previous one */
		Float angleThreshold = degToRad(props.getFloat("angleThreshold", 1.0f));
		Float dpThresh = std::cos(angleThreshold);

		/* Object-space -> World-space transformation */
		Transform objectToWorld = props.getTransform("toWorld", Transform());

		Log(EInfo, "Loading hair geometry from \"%s\" ..", path.leaf().c_str());

		fs::ifstream is(path);
		if (is.fail())
			Log(EError, "Could not open \"%s\"!", path.file_string().c_str());

		std::string line;
		bool newFiber = true;
		Point p, lastP(Point::Zero());
		std::vector<Point> vertices;
		std::vector<bool> vertexStartsFiber;
		Vector tangent(Vector::Zero());
		size_t nDegenerate = 0, nSkipped = 0;

		while (is.good()) {
			std::getline(is, line);
			if (line.norm() > 0 && line[0] == '#') {
				newFiber = true;
				continue;
			}
			std::istringstream iss(line);
			iss >> p.x >> p.y >> p.z;
			if (!iss.fail()) {
				p = objectToWorld(p);
				if (newFiber) {
					vertices.push_back(p);
					vertexStartsFiber.push_back(newFiber);
					lastP = p;
					tangent = Vector::Zero();
				} else if (p != lastP) {
					if (tangent.isZero()) {
						vertices.push_back(p);
						vertexStartsFiber.push_back(newFiber);
						tangent = (p - lastP).normalized();
						lastP = p;
					} else {
						Vector nextTangent = (p - lastP).normalized();
						if (nextTangent.dot(tangent) > dpThresh) {
							/* Too small of a difference in the tangent value,
							   just overwrite the previous vertex by the current one */
							tangent = (p - vertices[vertices.size()-2]).normalized();
							vertices[vertices.size()-1] = p;
							++nSkipped;
						} else {
							vertices.push_back(p);
							vertexStartsFiber.push_back(newFiber);
							tangent = nextTangent;
						}
						lastP = p;
					}
				} else {
					nDegenerate++;
				}
				newFiber = false;
			} else {
				newFiber = true;
			}
		}

		if (nDegenerate > 0)
			Log(EInfo, "Encountered " SIZE_T_FMT 
				" degenerate segments!", nDegenerate);
		if (nSkipped > 0)
			Log(EInfo, "Skipped " SIZE_T_FMT 
				" low-curvature segments.", nSkipped);

		vertexStartsFiber.push_back(true);

		m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
	}

	Hair(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		Float radius = stream->readFloat();
		size_t vertexCount = (size_t) stream->readUInt();

		std::vector<Point> vertices(vertexCount);
		std::vector<bool> vertexStartsFiber(vertexCount+1);
		stream->readFloatArray((Float *) &vertices[0], vertexCount * 3);

		for (size_t i=0; i<vertexCount; ++i) 
			vertexStartsFiber[i] = stream->readBool();
		vertexStartsFiber[vertexCount] = true;

		m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		const std::vector<Point> &vertices = m_kdtree->getVertices();
		const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();

		stream->writeFloat(m_kdtree->getRadius());
		stream->writeUInt((uint32_t) vertices.size());
		stream->writeFloatArray((Float *) &vertices[0], vertices.size() * 3);
		for (size_t i=0; i<vertices.size(); ++i)
			stream->writeBool(vertexStartsFiber[i]);
	}

	bool rayIntersect(const Ray &ray, Float mint, 
			Float maxt, Float &t, void *temp) const {
		return m_kdtree->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		return m_kdtree->rayIntersect(ray, mint, maxt);
	}

	void fillIntersectionRecord(const Ray &ray, 
		const void *temp, Intersection &its) const {
		its.p = ray(its.t);

		/* No UV coordinates for now */
		its.uv = Point2::Zero();
		its.dpdu = Vector::Zero();
		its.dpdv = Vector::Zero();

		const HairKDTree::IndexType *storage = 
			static_cast<const HairKDTree::IndexType *>(temp);
		HairKDTree::IndexType iv = *storage;

		const Vector axis = m_kdtree->tangent(iv);
		its.geoFrame.s = axis;
		const Vector relHitPoint = its.p - m_kdtree->firstVertex(iv);
		its.geoFrame.n = Normal((relHitPoint - axis.dot(relHitPoint) * axis).normalized());
		its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.hasUVPartials = false;
		its.shape = this;
	}

	ref<TriMesh> createTriMesh() {
		size_t nSegments = m_kdtree->getSegmentCount();
		/// Use very approximate geometry for large hair meshes
		const size_t phiSteps = (nSegments > 100000) ? 4 : 10;
		const Float dPhi   = (2*M_PI) / phiSteps;

		ref<TriMesh> mesh = new TriMesh("Hair mesh approximation",
			phiSteps*2*nSegments, phiSteps*2*nSegments, true, false, false);

		Point *vertices = mesh->getVertexPositions();
		Normal *normals = mesh->getVertexNormals();
		Triangle *triangles = mesh->getTriangles();
		size_t triangleIdx = 0, vertexIdx = 0;
		
		const std::vector<Point> &hairVertices = m_kdtree->getVertices();
		const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();
		const Float radius = m_kdtree->getRadius();
		Float *cosPhi = new Float[phiSteps];
		Float *sinPhi = new Float[phiSteps];
		for (size_t i=0; i<phiSteps; ++i) {
			sinPhi[i] = std::sin(i*dPhi);
			cosPhi[i] = std::cos(i*dPhi);
		}

		size_t hairIdx = 0;
		for (size_t iv=0; iv<hairVertices.size()-1; iv++) {
			if (!vertexStartsFiber[iv+1]) {
				for (size_t phi=0; phi<phiSteps; ++phi) {
					Vector tangent = m_kdtree->tangent(iv);
					Vector dir = Frame(tangent).toWorld(
							Vector(cosPhi[phi], sinPhi[phi], 0));
					Normal miterNormal1 = m_kdtree->firstMiterNormal(iv);
					Normal miterNormal2 = m_kdtree->secondMiterNormal(iv);
					Float t1 = radius * miterNormal1.dot(dir) / miterNormal1.dot(tangent);
					Float t2 = radius * miterNormal2.dot(dir) / miterNormal2.dot(tangent);

					Normal normal(dir.normalized());
					normals[vertexIdx] = normal;
					vertices[vertexIdx++] = m_kdtree->firstVertex(iv) + radius*dir - tangent*t1;
					normals[vertexIdx] = normal;
					vertices[vertexIdx++] = m_kdtree->secondVertex(iv) + radius*dir - tangent*t2;

					int idx0 = 2*(phi + hairIdx*phiSteps), idx1 = idx0+1;
					int idx2 = (2*phi+2) % (2*phiSteps) + 2*hairIdx*phiSteps, idx3 = idx2+1;
					triangles[triangleIdx].idx[0] = idx0;
					triangles[triangleIdx].idx[1] = idx2;
					triangles[triangleIdx].idx[2] = idx1;
					triangleIdx++;
					triangles[triangleIdx].idx[0] = idx1;
					triangles[triangleIdx].idx[1] = idx2;
					triangles[triangleIdx].idx[2] = idx3;
					triangleIdx++;
				}
				hairIdx++;
			}
		}
		Assert(triangleIdx == phiSteps*2*nSegments);
		Assert(vertexIdx == phiSteps*2*nSegments);

		delete[] cosPhi;
		delete[] sinPhi;

		mesh->setBSDF(m_bsdf);
		mesh->setLuminaire(m_luminaire);
		mesh->configure();

		return mesh.get();
	}

	const KDTreeBase<BoundingBox3> *getKDTree() const {
		return m_kdtree.get();
	}

	BoundingBox3 getBoundingBox() const {
		return m_kdtree->getBoundingBox();
	}

	Float getSurfaceArea() const {
		Log(EError, "Hair::getSurfaceArea(): Not implemented.");
		return -1;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Hair[" << endl
			<< "   numVertices = " << m_kdtree->getVertexCount() << ","
			<< "   numSegments = " << m_kdtree->getSegmentCount() << ","
			<< "   numHairs = " << m_kdtree->getHairCount() << ","
			<< "   radius = " << m_kdtree->getRadius()
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<HairKDTree> m_kdtree;
};

MTS_IMPLEMENT_CLASS(HairKDTree, false, KDTreeBase)
MTS_IMPLEMENT_CLASS_S(Hair, false, Shape)
MTS_EXPORT_PLUGIN(Hair, "Hair intersection primitive");
MTS_NAMESPACE_END
