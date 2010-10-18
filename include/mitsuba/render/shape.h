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

#if !defined(__SHAPE_H)
#define __SHAPE_H

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/** \brief Data record, which holds sampling-related information
 *  for a shape.
 */
struct MTS_EXPORT_RENDER ShapeSamplingRecord {
public:
	/// Create a sampling record (does no initialization!)
	inline ShapeSamplingRecord() { }

	/// Create a sampling record with the specified position and normal
	inline ShapeSamplingRecord(const Point &p, const Normal &n)
		: p(p), n(n) { }

	/// Return a string representation
	std::string toString() const;
public:
	/// Sampled surface position
	Point p;

	/// Sampled surface normal
	Normal n;
};

/** \brief Container for all information related to
 * a surface intersection
 */
struct MTS_EXPORT_RENDER Intersection {
public:
	inline Intersection() :
		t(std::numeric_limits<Float>::infinity()), 
		shape(NULL) { }

	/* Convert a vector expressed inside the shading frame into world
	   coordinates */
	inline Vector toWorld(const Vector &v) const {
		return shFrame.toWorld(v);
	}

	/* Convert a vector expressed inside world coordinates frame into 
	   shading frame coordinates */
	inline Vector toLocal(const Vector &v) const {
		return shFrame.toLocal(v);
	}

	/// Is the current intersection valid?
	inline bool isValid() const {
		return t != std::numeric_limits<Float>::infinity();
	}

	/// Is the intersected shape also a luminaire?
	inline bool isLuminaire() const;

	/// Does the intersected shape have a subsurface integrator?
	inline bool hasSubsurface() const;

	/**
	 * Returns the BSDF of the intersected shape. The
	 * parameter ray must match the one used to create
	 * the intersection record. Computes texture coordinate
	 * partials if this is required by the BSDF.
	 * Should only be called if there is a valid
	 * intersection!
	 */
	inline const BSDF *getBSDF(const RayDifferential &ray);

	/**
	 * Returns radiance emitted into direction d.
	 * Should only be called if the intersected
	 * shape is indeed a luminaire!
	 */
	inline Spectrum Le(const Vector &d) const;

	/**
	 * Returns radiance from a subsurface integrator
	 * emitted into direction d.
	 * Should only be called if the intersected
	 * shape does indeed have a subsurface integrator!
	 */
	inline Spectrum LoSub(const Scene *scene, const Vector &d) const;

	/// Computes texture coordinate partials
	void computePartials(const RayDifferential &ray);

	/// Return a string representation
	std::string toString() const;
public:
	/// Distance traveled along the ray
	Float t;

	/* Intersection point in 3D coordinates */
	Point p;

	/// Geometry frame
	Frame geoFrame;

	/// Shading frame
	Frame shFrame;

	/// UV surface coordinates
	Point2 uv;

	/// Position partials wrt. the texture space parameterization
	Vector dpdu, dpdv;

	/// Texture coordinate mapping partials wrt. changes in screen-space
	Float dudx, dudy, dvdx, dvdy;

	/// Interpolated vertex color
	Spectrum color;

	/// Incident direction in the local frame
	Vector wi;

	/// Affected shape
	const Shape *shape;

	/// Have texture coordinate partials been computed
	bool hasUVPartials;
};

/** \brief Abstract base class of all shapes
 */
class MTS_EXPORT_RENDER Shape : public ConfigurableObject {
public:
	/// Is this a compound shape consisting of several sub-objects?
	virtual bool isCompound() const;

	/**
	 * \brief Return a sub-element of a compound shape. 
	 *
	 * When expanding shapes, the scene will repeatedly call this
	 * function with increasing indices. Returning \a NULL indicates
	 * that no more elements are available.
	 */
	virtual Shape *getElement(int i);

	/// Return the shape's surface area
	virtual Float getSurfaceArea() const = 0;

	/// Return a bounding box containing the shape
	virtual AABB getAABB() const = 0;

	/**
	 * \brief Returns the minimal axis-aligned bounding box 
	 * of this shape when clipped to another bounding box.
	 * 
	 * This is extremely important to construct decent kd-trees.
	 * The default implementation just takes the bounding box
	 * returned by \ref getAABB() and clips it to \a box.
	 */
	virtual AABB getClippedAABB(const AABB &box) const;

	/**
	 * \brief Fast ray intersection test
	 *
	 * Check whether the shape is intersected by the given ray. Some
	 * temporary space (\ref MTS_KD_INTERSECTION_TEMP-4 bytes) is,
	 * supplied which can be used to cache information about the 
	 * intersection. The function \ref fillIntersectionRecord() 
	 * can later use this information to fill in a detailed 
	 * intersection record.
	 */
	virtual bool rayIntersect(const Ray &ray, Float mint, 
			Float maxt, Float &t, void *temp) const;

	/**
	 * \brief Given that an intersection has been found, create a 
	 * detailed intersection record
	 */
	virtual void fillIntersectionRecord(const Ray &ray, Float t, 
			const void *temp, Intersection &its) const;

	/**
	 * \brief Sample a point on the shape
	 *
	 * Should be uniform wrt. surface area. Returns the 
	 * associated probability density
	 */
	virtual Float sampleArea(ShapeSamplingRecord &sRec, 
			const Point2 &sample) const;

	/**
	 * Return the probability density of sampling the 
	 * given point using sampleArea()
	 */
	virtual Float pdfArea(const ShapeSamplingRecord &sRec) const;

	/**
	 * \brief Sample a point on the shape and return the associated
	 * probability wrt. solid angle.
	 *
	 * Should ideally be uniform wrt. solid angle as seen 
	 * from \a x. The default implementation, just uses
	 * \ref sampleArea, which can lead to lots of noise.
	 */
	virtual Float sampleSolidAngle(ShapeSamplingRecord &sRec, 
			const Point &x, const Point2 &sample) const;

	/**
	 * \brief Return the probability density of sampling the given 
	 * point using \ref sampleSolidAngle().
	 */
	virtual Float pdfSolidAngle(const ShapeSamplingRecord &sRec, 
			const Point &x) const;

	/// Return the shape's BSDF
	inline const BSDF *getBSDF() const { return m_bsdf.get(); }
	/// Return the shape's BSDF
	inline BSDF *getBSDF() { return m_bsdf.get(); }

	/// Return the name of this shape
	virtual std::string getName() const;

	/// Does this shape have a sub-surface integrator?
	inline bool hasSubsurface() const { return m_subsurface.get() != NULL; }
	/// Return the associated sub-surface integrator
	inline Subsurface *getSubsurface() { return m_subsurface; }
	/// Return the associated sub-surface integrator 
	inline const Subsurface *getSubsurface() const { return m_subsurface.get(); }

	/// Is this shape also an area luminaire?
	inline bool isLuminaire() const { return m_luminaire.get() != NULL; }
	/// Return the associated luminaire (if any)
	inline Luminaire *getLuminaire() { return m_luminaire; }
	/// Return the associated luminaire (if any)
	inline const Luminaire *getLuminaire() const { return m_luminaire.get(); }

	/// Called once after parsing
	virtual void configure();
	/// Serialize this shape to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Add a child (e.g. a luminaire/sub surface integrator) to this shape
	void addChild(const std::string &name, ConfigurableObject *child);

	MTS_DECLARE_CLASS()
protected:
	/// Create a new shape
	Shape(const Properties &props);

	/// Unserialize a shape
	Shape(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Shape();
protected:
	ref<BSDF> m_bsdf;
	ref<Subsurface> m_subsurface;
	ref<Luminaire> m_luminaire;
};

MTS_NAMESPACE_END

#endif /* __SHAPE_H */


