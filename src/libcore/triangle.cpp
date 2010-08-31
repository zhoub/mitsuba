#include <mitsuba/core/triangle.h>

MTS_NAMESPACE_BEGIN

bool Triangle::rayIntersect(const Vertex *buffer, const Ray &ray,
	Float &u, Float &v, Float &t) const {
	const Point &v0 = buffer[idx[0]].v;
	const Point &v1 = buffer[idx[1]].v;
	const Point &v2 = buffer[idx[2]].v;

	/* find vectors for two edges sharing v[0] */
	Vector edge1 = v1 - v0, edge2 = v2 - v0;

	/* begin calculating determinant - also used to calculate U parameter */
	Vector pvec = cross(ray.d, edge2);

	/* if determinant is near zero, ray lies in plane of triangle */
	Float det = dot(edge1, pvec);

	if (det > -Epsilon && det < Epsilon)
		return false;
	Float inv_det = 1.0f / det;

	/* calculate distance from v[0] to ray origin */
	Vector tvec = ray.o - v0;

	/* calculate U parameter and test bounds */
	u = dot(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	/* prepare to test V parameter */
	Vector qvec = cross(tvec, edge1);

	/* calculate V parameter and test bounds */
	v = dot(ray.d, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	/* calculate t, ray intersects triangle */
	t = dot(edge2, qvec) * inv_det;

	return true;
}

Point Triangle::sample(const Vertex *buffer, Normal &normal, 
	const Point2 &sample) const {
	const Point &v0 = buffer[idx[0]].v;
	const Point &v1 = buffer[idx[1]].v;
	const Point &v2 = buffer[idx[2]].v;
	const Normal &n0 = buffer[idx[0]].n;
	const Normal &n1 = buffer[idx[1]].n;
	const Normal &n2 = buffer[idx[2]].n;

	Point2 bary = squareToTriangle(sample);
	Vector sideA = v1 - v0, sideB = v2 - v0;
	Point p = v0 + (sideA * bary.x) + (sideB * bary.y);	
	normal = Normal(normalize(
		n0 * (1.0f - bary.x - bary.y) +
		n1 * bary.x + n2 * bary.y
	));

	return p;
}

Float Triangle::surfaceArea(const Vertex *buffer) const {
	const Point &v0 = buffer[idx[0]].v;
	const Point &v1 = buffer[idx[1]].v;
	const Point &v2 = buffer[idx[2]].v;
	Vector sideA = v1 - v0, sideB = v2 - v0;
	return 0.5f * cross(sideA, sideB).length();
}

AABB Triangle::getAABB(const Vertex *buffer) const {
	AABB aabb;
	for (int k=0; k<3; k++)
		aabb.expandBy(buffer[idx[k]].v);
	return aabb;
}

inline void sutherlandHodgman(std::vector<Point> &vertices, int axis, 
	Float splitPos, bool isMinimum) {
	int vertexCount = (int) vertices.size();
	if (vertexCount < 2)
		return;

	Point cur = vertices[0];
	Float sign = isMinimum ? 1.0f : -1.0f, 
	      distance = sign * (cur[axis] - splitPos);
	bool curIsInside = (distance >= 0);

	for (int i=0; i<vertexCount; i++) {
		Point next = vertices[(i+1)%vertexCount];
		distance = sign * (next[axis] - splitPos);
		bool nextIsInside = (distance >= 0);

		if (curIsInside && nextIsInside) {
			/* Both this and the next vertex are inside, add to the list */
			vertices.push_back(next);
		} else if (curIsInside && !nextIsInside) {
			/* Going outside -- add the intersection */
			Float t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			vertices.push_back(cur + (next - cur) * t);
		} else if (!curIsInside && nextIsInside) {
			/* Coming back inside -- add the intersection + next vertex */
			Float t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			vertices.push_back(cur + (next - cur) * t);
			vertices.push_back(next);
		} else {
			/* Entirely outside - do not add anything */
		}
		cur = next;
		curIsInside = nextIsInside;
	}
	vertices.erase(vertices.begin(), vertices.begin() + vertexCount);
}

AABB Triangle::getClippedAABB(const Vertex *buffer, const AABB &aabb) const {
	std::vector<Point> vertices;
	/* Reserve room for some additional vertices */
	vertices.reserve(8);
	for (int i=0; i<3; ++i)
		vertices.push_back(buffer[idx[i]].v);

	for (int axis=0; axis<3; ++axis) {
		sutherlandHodgman(vertices, axis, aabb.min[axis], true);
		sutherlandHodgman(vertices, axis, aabb.max[axis], false);
	}

	AABB result;
	for (unsigned int i=0; i<vertices.size(); ++i)
		result.expandBy(vertices[i]);

	/* Cover up some numerical imprecisions */

	for (int i=0; i<3; ++i)
		result.min[i] -= Epsilon * std::abs(result.min[i]);
	for (int i=0; i<3; ++i)
		result.max[i] += Epsilon * std::abs(result.max[i]);

	result.clip(aabb);

	return result;
}

MTS_NAMESPACE_END
