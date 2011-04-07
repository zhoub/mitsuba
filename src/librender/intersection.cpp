#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

void Intersection::computePartials(const RayDifferential &ray) {
	Float A[2][2], Bx[2], By[2], x[2];
	int axes[2];

	/* Compute the texture coordinates partials wrt. 
	   changes in the screen-space position. Based on PBRT */
	if (hasUVPartials)
		return;
	hasUVPartials = true;

	if (!ray.hasDifferentials || (dpdu.isZero() && dpdv.isZero())) {
		dudx = dvdx = dudy = dvdy = 0.0f;
		return;
	}

	/* Offset of the plane passing through the surface */
	const Float d = -dot(geoFrame.n, Vector(p));

	const Float txRecip = dot(geoFrame.n, ray.rx.d),
				tyRecip = dot(geoFrame.n, ray.ry.d);

	if (EXPECT_NOT_TAKEN(txRecip == 0 || tyRecip == 0)) {
		dudx = dvdx = dudy = dvdy = 0.0f;
		return;
	}

	/* Ray distances traveled */
	const Float tx = -(dot(geoFrame.n, Vector(ray.rx.o)) + d) / 
		txRecip;
	const Float ty = -(dot(geoFrame.n, Vector(ray.ry.o)) + d) / 
		tyRecip;

	/* Calculate the U and V partials by solving two out
	   of a set of 3 equations in an overconstrained system */
	Float absX = std::abs(geoFrame.n.x),
		  absY = std::abs(geoFrame.n.y),
		  absZ = std::abs(geoFrame.n.z);

	if (absX > absY && absX > absZ) {
		axes[0] = 1; axes[1] = 2;
	} else if (absY > absZ) {
		axes[0] = 0; axes[1] = 2;
	} else {
		axes[0] = 0; axes[1] = 1;
	}

	A[0][0] = dpdu[axes[0]];
	A[0][1] = dpdv[axes[0]];
	A[1][0] = dpdu[axes[1]];
	A[1][1] = dpdv[axes[1]];

	/* Auxilary intersection point of the adjacent rays */
	Point px = ray.rx(tx), py = ray.ry(ty);
	Bx[0] = px[axes[0]] - p[axes[0]];
	Bx[1] = px[axes[1]] - p[axes[1]];
	By[0] = py[axes[0]] - p[axes[0]];
	By[1] = py[axes[1]] - p[axes[1]];
	
	if (EXPECT_TAKEN(solveLinearSystem2x2(A, Bx, x))) {
		dudx = x[0]; dvdx = x[1];
	} else {
		dudx = 1; dvdx = 0;
	}

	if (EXPECT_TAKEN(solveLinearSystem2x2(A, By, x))) {
		dudy = x[0]; dvdy = x[1];
	} else {
		dudy = 0; dudy = 1;
	}
}

std::string Intersection::toString() const {
	if (!isValid())
		return "Intersection[invalid]";
	std::ostringstream oss;
	oss << "Intersection[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  t = " << t << "," << std::endl
		<< "  geoFrame = " << indent(geoFrame.toString()) << "," << std::endl
		<< "  shFrame = " << indent(shFrame.toString()) << "," << std::endl
		<< "  uv = " << uv.toString() << "," << std::endl
		<< "  hasUVPartials = " << hasUVPartials << "," << std::endl
		<< "  dpdu = " << dpdu.toString() << "," << std::endl
		<< "  dpdv = " << dpdv.toString() << "," << std::endl
		<< "  time = " << time << "," << std::endl
		<< "  shape = " << indent(((Object *)shape)->toString()) << std::endl
		<< "]";
	return oss.str();
}

MTS_NAMESPACE_END
