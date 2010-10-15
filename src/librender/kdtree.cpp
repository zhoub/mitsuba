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

#include <mitsuba/render/kdtree.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

KDTree::KDTree() {
#if !defined(MTS_KD_CONSERVE_MEMORY)
	m_triAccel = NULL;
#endif
	m_shapeMap.push_back(0);
}

KDTree::~KDTree() {
#if !defined(MTS_KD_CONSERVE_MEMORY)
	if (m_triAccel)
		freeAligned(m_triAccel);
#endif
	for (size_t i=0; i<m_shapes.size(); ++i)
		m_shapes[i]->decRef();
}

static StatsCounter raysTraced("General", "Normal rays traced");
static StatsCounter shadowRaysTraced("General", "Shadow rays traced");

void KDTree::addShape(const Shape *shape) {
	Assert(!isBuilt());
	if (shape->isCompound())
		Log(EError, "Cannot add compound shapes to a kd-tree - expand them first!");
	if (shape->getClass()->derivesFrom(TriMesh::m_theClass)) {
		// Triangle meshes are expanded into individual primitives,
		// which are visible to the tree construction code. Generic
		// primitives are only handled by their AABBs
		m_shapeMap.push_back((size_type) 
			static_cast<const TriMesh *>(shape)->getTriangleCount());
		m_triangleFlag.push_back(true);
	} else {
		m_shapeMap.push_back(1);
		m_triangleFlag.push_back(false);
	}
	shape->incRef();
	m_shapes.push_back(shape);
}

void KDTree::build() {
	for (size_t i=1; i<m_shapeMap.size(); ++i)
		m_shapeMap[i] += m_shapeMap[i-1];

	buildInternal();

#if !defined(MTS_KD_CONSERVE_MEMORY)
	ref<Timer> timer = new Timer();
	size_type primCount = getPrimitiveCount();
	Log(EDebug, "Precomputing triangle intersection information (%s)",
			memString(sizeof(TriAccel)*primCount).c_str());
	m_triAccel = static_cast<TriAccel *>(allocAligned(primCount * sizeof(TriAccel)));

	index_type idx = 0;
	for (index_type i=0; i<m_shapes.size(); ++i) {
		const Shape *shape = m_shapes[i];
		if (m_triangleFlag[i]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			const Triangle *triangles = mesh->getTriangles();
			const Point *positions = mesh->getVertexPositions();
			for (index_type j=0; j<mesh->getTriangleCount(); ++j) {
				const Triangle &tri = triangles[j];
				const Point &v0 = positions[tri.idx[0]];
				const Point &v1 = positions[tri.idx[1]];
				const Point &v2 = positions[tri.idx[2]];
				m_triAccel[idx].load(v0, v1, v2);
				m_triAccel[idx].shapeIndex = i;
				m_triAccel[idx].index = j;
				++idx;
			}
		} else {
			m_triAccel[idx].shapeIndex = i;
			m_triAccel[idx].k = KNoTriangleFlag;
			++idx;
		}
	}
	Log(EDebug, "Finished -- took %i ms.", timer->getMilliseconds());
	Log(EDebug, "");
	KDAssert(idx == primCount);
#endif
}

bool KDTree::rayIntersect(const Ray &ray, Intersection &its) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	its.t = std::numeric_limits<Float>::infinity(); 
	Float mint, maxt;

	++raysTraced;
	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		/* Use an adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::abs(ray.o.x), 
				std::abs(ray.o.y)), std::abs(ray.o.z));

		if (rayMinT > mint) mint = rayMinT;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) {

			if (rayIntersectHavran<false>(ray, mint, maxt, its.t, temp)) {
				/* After having found a unique intersection, fill a proper record
				   using the temporary information collected in \ref intersect() */
				const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(temp);
				const TriMesh *trimesh = static_cast<const TriMesh *>(m_shapes[cache->shapeIndex]);
				const Triangle &tri = trimesh->getTriangles()[cache->index];
				const Point *vertexPositions = trimesh->getVertexPositions();
				const Normal *vertexNormals = trimesh->getVertexNormals();
				const Point2 *vertexTexcoords = trimesh->getVertexTexcoords();
				const Spectrum *vertexColors = trimesh->getVertexColors();
				const TangentSpace *vertexTangents = trimesh->getVertexTangents();
				const Vector b(1 - cache->u - cache->v, cache->u, cache->v);

				const uint32_t idx0 = tri.idx[0], idx1 = tri.idx[1], idx2 = tri.idx[2];
				const Point &p0 = vertexPositions[idx0];
				const Point &p1 = vertexPositions[idx1];
				const Point &p2 = vertexPositions[idx2];

				//its.p = ray(its.t);
				its.p = p0 * b.x + p1 * b.y + p2 * b.z;
				Normal faceNormal(cross(p1-p0, p2-p0));
				Float length = faceNormal.length();
				if (!faceNormal.isZero())
					faceNormal /= length;

				its.geoFrame = Frame(faceNormal);

				if (EXPECT_TAKEN(vertexNormals)) {
					const Normal &n0 = vertexNormals[idx0];
					const Normal &n1 = vertexNormals[idx1];
					const Normal &n2 = vertexNormals[idx2];

					if (EXPECT_TAKEN(!vertexTangents)) {
						its.shFrame = Frame(normalize(n0 * b.x + n1 * b.y + n2 * b.z));
					} else {
						const TangentSpace &t0 = vertexTangents[idx0];
						const TangentSpace &t1 = vertexTangents[idx1];
						const TangentSpace &t2 = vertexTangents[idx2];
						const Vector dpdu = t0.dpdu * b.x + t1.dpdu * b.y + t2.dpdu * b.z;
						its.shFrame.n = normalize(n0 * b.x + n1 * b.y + n2 * b.z);
						its.shFrame.s = normalize(dpdu - its.shFrame.n 
							* dot(its.shFrame.n, dpdu));
						its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
						its.dpdu = dpdu;
						its.dpdv = t0.dpdv * b.x + t1.dpdv * b.y + t2.dpdv * b.z;
					}
				} else {
					its.shFrame = its.geoFrame;
				}

				if (EXPECT_TAKEN(vertexTexcoords)) {
					const Point2 &t0 = vertexTexcoords[idx0];
					const Point2 &t1 = vertexTexcoords[idx0];
					const Point2 &t2 = vertexTexcoords[idx0];
					its.uv = t0 * b.x + t1 * b.y + t2 * b.z;
				} else {
					its.uv = Point2(0.0f);
				}

				if (EXPECT_NOT_TAKEN(vertexColors)) {
					const Spectrum &c0 = vertexColors[idx0],
								   &c1 = vertexColors[idx1],
								   &c2 = vertexColors[idx2];
					its.color = c0 * b.x + c1 * b.y + c2 * b.z;
				}

				its.wi = its.toLocal(-ray.d);
				its.shape = trimesh;
				its.hasUVPartials = false;
				return true;
			}
		}
	}
	return false;
}

bool KDTree::rayIntersect(const Ray &ray) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	Float mint, maxt, t = std::numeric_limits<Float>::infinity();

	++shadowRaysTraced;
	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		/* Use an adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::abs(ray.o.x), 
				std::abs(ray.o.y)), std::abs(ray.o.z));

		if (rayMinT > mint) mint = rayMinT;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) 
			if (rayIntersectHavran<true>(ray, mint, maxt, t, temp)) 
				return true;
	}
	return false;
}


#if defined(MTS_HAS_COHERENT_RT)
static StatsCounter coherentPackets("General", "Coherent ray packets");
static StatsCounter incoherentPackets("General", "Incoherent ray packets");

void KDTree::rayIntersectPacket(const RayPacket4 &packet, 
		const RayInterval4 &rayInterval, Intersection4 &its) const {
	CoherentKDStackEntry MM_ALIGN16 stack[MTS_KD_MAXDEPTH];
	RayInterval4 MM_ALIGN16 interval;

	const KDNode * __restrict currNode = m_nodes;
	int stackIndex = 0;

	++coherentPackets;

	/* First, intersect with the kd-tree AABB to determine
	   the intersection search intervals */
	if (!m_aabb.rayIntersectPacket(packet, interval))
		return;

	interval.mint.ps = _mm_max_ps(interval.mint.ps, rayInterval.mint.ps);
	interval.maxt.ps = _mm_min_ps(interval.maxt.ps, rayInterval.maxt.ps);

	__m128 itsFound = _mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps),
		   masked = itsFound;
	if (_mm_movemask_ps(itsFound) == 0xF)
		return;

	while (currNode != NULL) {
		while (EXPECT_TAKEN(!currNode->isLeaf())) {
			const uint8_t axis = currNode->getAxis();

			/* Calculate the plane intersection */
			const __m128
				splitVal = _mm_set1_ps(currNode->getSplit()),
				t = _mm_mul_ps(_mm_sub_ps(splitVal, packet.o[axis].ps),
					packet.dRcp[axis].ps);

			const __m128
				startsAfterSplit = _mm_or_ps(masked, 
					_mm_cmplt_ps(t, interval.mint.ps)),
				endsBeforeSplit = _mm_or_ps(masked,
					_mm_cmpgt_ps(t, interval.maxt.ps));

			currNode = currNode->getLeft() + packet.signs[axis][0];

			/* The interval completely completely lies on one side
			   of the split plane */
			if (EXPECT_TAKEN(_mm_movemask_ps(startsAfterSplit) == 15)) {
				currNode = currNode->getSibling();
				continue;
			}

			if (EXPECT_TAKEN(_mm_movemask_ps(endsBeforeSplit) == 15)) 
				continue;

			stack[stackIndex].node = currNode->getSibling();
			stack[stackIndex].interval.maxt =    interval.maxt;
			stack[stackIndex].interval.mint.ps = _mm_max_ps(t, interval.mint.ps);
			interval.maxt.ps =                   _mm_min_ps(t, interval.maxt.ps);
			masked = _mm_or_ps(masked, 	
					_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
			stackIndex++;
		}

		/* Arrived at a leaf node - intersect against primitives */
		const index_type primStart = currNode->getPrimStart();
		const index_type primEnd = currNode->getPrimEnd();

		if (EXPECT_NOT_TAKEN(primStart != primEnd)) {
#ifdef MTS_USE_TRIACCEL4
			const int count = m_packedTriangles[primStart].indirectionCount;
			primStart = m_packedTriangles[primStart].indirectionIndex;
			primEnd = primStart + count;
#endif
			__m128 
				searchStart = _mm_max_ps(rayInterval.mint.ps, 
					_mm_mul_ps(interval.mint.ps, SSEConstants::om_eps.ps)),
				searchEnd   = _mm_min_ps(rayInterval.maxt.ps, 
					_mm_mul_ps(interval.maxt.ps, SSEConstants::op_eps.ps));

			for (index_type entry=primStart; entry != primEnd; entry++) {
				const TriAccel &kdTri = m_triAccel[m_indices[entry]];
				if (EXPECT_TAKEN(kdTri.k != KNoTriangleFlag)) {
					itsFound = _mm_or_ps(itsFound, 
						kdTri.rayIntersectPacket(packet, searchStart, searchEnd, masked, its));
				} else {
					/* Not a triangle - invoke the shape's intersection routine */
					__m128 hasIts = m_shapes[kdTri.shapeIndex]->rayIntersectPacket(packet, 
							searchStart, searchEnd, masked, its);
					itsFound = _mm_or_ps(itsFound, hasIts);
					its.primIndex.pi  = mux_epi32(pstoepi32(hasIts), 
						load1_epi32(kdTri.index), its.primIndex.pi);
					its.shapeIndex.pi = mux_epi32(pstoepi32(hasIts),
						load1_epi32(kdTri.shapeIndex), its.shapeIndex.pi);
				}
				searchEnd = _mm_min_ps(searchEnd, its.t.ps);
			}
		}

		/* Abort if the tree has been traversed or if
		   intersections have been found for all four rays */
		if (_mm_movemask_ps(itsFound) == 0xF || --stackIndex < 0)
			break;

		/* Pop from the stack */
		currNode = stack[stackIndex].node;
		interval = stack[stackIndex].interval;
		masked = _mm_or_ps(itsFound, 
			_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
	}
}
	
void KDTree::rayIntersectPacket(const Ray *rays, Intersection *itsArray) const {
	RayPacket4 MM_ALIGN16 packet;
	RayInterval4 MM_ALIGN16 interval(rays);
	Intersection4 MM_ALIGN16 its4;

	if (packet.load(rays)) {
		rayIntersectPacket(packet, interval, its4);

		for (int i=0; i<4; i++) {
			Intersection &its = itsArray[i];

			its.t = its4.t.f[i];

			if (its.t != std::numeric_limits<float>::infinity()) {
				const uint32_t shapeIndex = its4.shapeIndex.i[i];
				const uint32_t primIndex = its4.primIndex.i[i];
				const Shape *shape = m_shapes[shapeIndex];

				if (EXPECT_TAKEN(primIndex != KNoTriangleFlag)) {
					const TriMesh *trimesh = static_cast<const TriMesh *>(shape);
					const Triangle &tri = trimesh->getTriangles()[primIndex];
					const Point *vertexPositions = trimesh->getVertexPositions();
					const Normal *vertexNormals = trimesh->getVertexNormals();
					const Point2 *vertexTexcoords = trimesh->getVertexTexcoords();
					const TangentSpace *vertexTangents = trimesh->getVertexTangents();
					const Vector b(1 - its4.u.f[i] - its4.v.f[i], its4.u.f[i], its4.v.f[i]);

					const uint32_t idx0 = tri.idx[0], idx1 = tri.idx[1], idx2 = tri.idx[2];
					const Point &p0 = vertexPositions[idx0];
					const Point &p1 = vertexPositions[idx1];
					const Point &p2 = vertexPositions[idx2];

					//its.p = ray(its.t);
					its.p = p0 * b.x + p1 * b.y + p2 * b.z;
					its.geoFrame = Frame(normalize(cross(p1-p0, p2-p0)));
					if (EXPECT_TAKEN(vertexNormals)) {
						const Normal &n0 = vertexNormals[idx0];
						const Normal &n1 = vertexNormals[idx1];
						const Normal &n2 = vertexNormals[idx2];

						if (EXPECT_TAKEN(!vertexTangents)) {
							its.shFrame = Frame(normalize(n0 * b.x + n1 * b.y + n2 * b.z));
						} else {
							const TangentSpace &t0 = vertexTangents[idx0];
							const TangentSpace &t1 = vertexTangents[idx1];
							const TangentSpace &t2 = vertexTangents[idx2];
							const Vector dpdu = t0.dpdu * b.x + t1.dpdu * b.y + t2.dpdu * b.z;
							its.shFrame.n = normalize(n0 * b.x + n1 * b.y + n2 * b.z);
							its.shFrame.s = normalize(dpdu - its.shFrame.n 
								* dot(its.shFrame.n, dpdu));
							its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
							its.dpdu = dpdu;
							its.dpdv = t0.dpdv * b.x + t1.dpdv * b.y + t2.dpdv * b.z;
						}
					} else {
						its.shFrame = its.geoFrame;
					}

					if (vertexTexcoords) {
						const Point2 &t0 = vertexTexcoords[idx0];
						const Point2 &t1 = vertexTexcoords[idx0];
						const Point2 &t2 = vertexTexcoords[idx0];
						its.uv = t0 * b.x + t1 * b.y + t2 * b.z;
					}

					its.wi = its.toLocal(-rays[i].d);
					its.hasUVPartials = false;
					its.shape = trimesh;
				} else {
					/* Non-triangle shape: intersect again to fill in details */
					shape->rayIntersect(rays[i], its);
				}
			}
		}
	} else {
		rayIntersectPacketIncoherent(rays, itsArray);
	}
}

void KDTree::rayIntersectPacketIncoherent(const RayPacket4 &packet, 
		const RayInterval4 &rayInterval, Intersection4 &its4) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];

	++incoherentPackets;
	for (int i=0; i<4; i++) {
		Ray ray;
		for (int axis=0; axis<3; axis++) {
			ray.o[axis] = packet.o[axis].f[i];
			ray.d[axis] = packet.d[axis].f[i];
			ray.dRcp[axis] = packet.dRcp[axis].f[i];
		}
		ray.mint = rayInterval.mint.f[i];
		ray.maxt = rayInterval.maxt.f[i];
		if (ray.mint < ray.maxt && rayIntersectHavran<false>(ray, ray.mint, ray.maxt, its4.t.f[i], temp)) {
			const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(temp);
			its4.u.f[i] = cache->u;
			its4.v.f[i] = cache->u;
			its4.shapeIndex.i[i] = cache->shapeIndex;
			its4.primIndex.i[i] = cache->index;
		}
	}
}

#else // !MTS_HAS_COHERENT_RT
void KDTree::rayIntersectPacket(const Ray *rays, Intersection *its) const {
	rayIntersectPacketIncoherent(rays, its);
}
#endif

MTS_IMPLEMENT_CLASS(KDTree, false, GenericKDTree)
MTS_NAMESPACE_END
