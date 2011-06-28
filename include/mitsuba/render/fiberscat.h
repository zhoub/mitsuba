/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__HAIRSCAT_H)
#define __HAIRSCAT_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/render/common.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Data structure, which contains information 
 * required to sample or query a hair scattering function. 
 */
struct MTS_EXPORT_RENDER FiberScatteringRecord {
	/**
	 * \brief Normalized incident direction vector, which points away
	 * from the scattering event.
	 */
	Vector wi;

	/// Normalized outgoing direction vector
	Vector wo;

	/* Transported quantity (radiance or importance) -- required for 
	   rendering with non-reciprocal scattering functions */
	ETransportQuantity quantity;

	inline FiberScatteringRecord(const Vector &wi) : 
		wi(wi), quantity(ERadiance) {
	}

	inline FiberScatteringRecord(const Vector &wi, 
			const Vector &wo) : wi(wi),  wo(wo), quantity(ERadiance) {
	}

	std::string toString() const;
};

/** \brief Abstract hair scattering function.
 */
class MTS_EXPORT_RENDER FiberScatteringFunction : public ConfigurableObject {
public:
	/**
	 * \brief Evaluate the hair scattering function for a 
	 * direction pair (wi, wo)
	 */
	virtual Spectrum f(const FiberScatteringRecord &hRec) const = 0;

	/**
	 * \brief Importance sample the hair scattering function. 
	 *
	 * \param sampler
	 *     Sample generator
	 * \return
	 *     Weight value equal to the throughput divided by 
	 *     the probability of the sampled direction.
	 */
	virtual Spectrum sample(FiberScatteringRecord &hRec, 
		Sampler *sampler) const = 0;

	/**
	 * \brief Compute the probability of sampling wo (given wi).
	 */
	virtual Float pdf(const FiberScatteringRecord &hRec) const = 0;

	/// Return a string representation
	virtual std::string toString() const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Create a new hair scattering function instance
	inline FiberScatteringFunction(const Properties &props) :
		ConfigurableObject(props) { }

	/// Unserialize a hair scattering function
	inline FiberScatteringFunction(Stream *stream, InstanceManager *manager) :
		ConfigurableObject(stream, manager) { }

	/// Virtual destructor
	virtual ~FiberScatteringFunction() { }
};

MTS_NAMESPACE_END

#endif /* __HAIRSCAT_H */
