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

#if !defined(__MEDIUM_H)
#define __MEDIUM_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN
/**
 * \brief Data record associated with the sampling procedure responsible for
 * choosing a point on the in-scattering line integral of the RTE
 *
 * \sa Medium::sampleDistance()
 */
struct MTS_EXPORT_RENDER MediumSamplingRecord {
public:
	inline MediumSamplingRecord() { }

	/// Return a string representation
	std::string toString() const;
public:

	/// Traveled distance
	Float t;

	/// Location of the scattering interaction
	Point p;

	/// Local particle orientation at \ref p
	Vector orientation;

	/**
	 * \brief Specifies the transmittance along the segment [mint, t]
	 *
	 * When sampling a distance fails, this contains the 
	 * transmittance along the whole ray segment [mint, maxDist].
	 */
	Spectrum transmittance;

	/// The medium's absorption coefficient at \ref p
	Spectrum sigmaA;

	/// The medium's scattering coefficient at \ref p
	Spectrum sigmaS;

	/// Records the probability of sampling a medium interaction at p
	Float pdfSuccess;
	
	/**
	 * \brief Records the probability of sampling a medium 
	 * interaction in the reverse direction
	 *
	 * This is essentially the density of obtained by calling \ref sampleDistance,
	 * but starting at \c p and stopping at \c ray.o. These probabilities
	 * are important for bidirectional methods.
	 */
	Float pdfSuccessRev;

	/**
	 * When the \ref Medium::sampleDistance() is successful, this function
	 * returns the probability of not having generated a medium interaction
	 * until \ref t. Otherwise, it records the probability of
	 * not generating any interactions in the whole interval [mint, maxt].
	 * This probability is assumed to be symmetric with respect to
	 * sampling from the other direction, which is why there is no
	 * \a pdfFailureRev field.
	 */
	Float pdfFailure;

	/// Max. single scattering albedo over all spectral samples
	Float albedo;
};

/** \brief Abstract participating medium 
 */
class MTS_EXPORT_RENDER Medium : public NetworkedObject {
public:
	/** 
	 * \brief Compute the transmittance along a ray segment
	 *
	 * Computes the transmittance along a ray segment 
	 * [mint, maxt] associated with the ray. It is assumed
	 * that the ray has a normalized direction value.
	 *
	 */
	virtual Spectrum getTransmittance(const Ray &ray) const = 0;

	/**
	 * \brief Sample a distance along the ray segment [mint, maxt]
	 *
	 * Should ideally importance sample with respect to the transmittance.
	 * It is assumed that the ray has a normalized direction value.
	 *
	 * \param ray      Ray, along which a distance should be sampled
	 * \param mRec     Medium sampling record to be filled with the result
	 * \return         \a false if the maximum distance was exceeded, or if
	 *                 no interaction inside the medium could be sampled.
	 */
	virtual bool sampleDistance(const Ray &ray,
		MediumSamplingRecord &mRec, Sampler *sampler) const = 0;

	/**
	 * \brief Compute the density of sampling distance \a t along the 
	 * ray using the sampling strategy implemented by \a sampleDistance. 
	 *
	 * The function computes the continuous densities in the case of
	 * a successful \ref sampleDistance() invocation (in both directions),
	 * as well as the Dirac delta density associated with a failure.
	 * For convenience, it also stores the transmittance along the ray
	 * segment in \a mRec.
	 */
	virtual void pdfDistance(const Ray &ray, Float t, 
		MediumSamplingRecord &mRec) const = 0;

	/// Return the phase function of this medium
	inline const PhaseFunction *getPhaseFunction() const { return m_phaseFunction.get(); }

	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	/// Serialize this medium to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Add a child ConfigurableObject
	virtual void addChild(const std::string &name, ConfigurableObject *child);

	/// Return a bounding volume
	inline const AABB &getAABB() const { return m_aabb; }

	/// Return a string representation
	virtual std::string toString() const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Create a new participating medium instance
	Medium(const Properties &props);
	
	/// Unserialize a participating medium
	Medium(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Medium() { }
protected:	
	Spectrum m_sigmaA;
	Spectrum m_sigmaS;
	Spectrum m_sigmaT;
	Float m_albedo;
	AABB m_aabb;
	ref<PhaseFunction> m_phaseFunction;
	Float m_densityMultiplier;
};

MTS_NAMESPACE_END

#endif /* __MEDIUM_H */
