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

#if !defined(__REPLAYABLE_SAMPLER_H)
#define __REPLAYABLE_SAMPLER_H

#include <mitsuba/render/sampler.h>

MTS_NAMESPACE_BEGIN

/**
 * Specialized sampler implementation used to seed MLT-style algorithm. Allows
 * to query for the current sample index, which can later be used to rewind
 * back to this state. In the case of MLT, this makes it possible to sample paths 
 * approximately proportional to their contribution without actually having
 * to store millions of path. Note that `rewinding' is naive -- it just 
 * resets & regenerates the whole random number sequence, which might be slow.
 */
class MTS_EXPORT_BIDIR ReplayableSampler : public Sampler {
public:
	/// Construct a new sampler
	ReplayableSampler();

	/// Unserialize a sampler
	ReplayableSampler(Stream *stream, InstanceManager *manager);

	/**
	 * Create a clone of this sampler. The clone is allowed to be different
	 * to some extent, e.g. a pseudorandom generator should be based on a
	 * different random seed compared to the original. All other parameters,
	 * are copied exactly.
	 */
	virtual ref<Sampler> clone();

	/* Does nothing in this implementation */
	virtual void advance();
	virtual void generate();

	/// Manually set the current sample index
	virtual void setSampleIndex(size_t sampleIndex);

	/// Retrieve the next component value from the current sample
	virtual Float next1D();

	/// Retrieve the next two component values from the current sample
	virtual Point2 next2D();

	/* Unsupported by this implementation */
	virtual void request2DArray(unsigned int size);
	virtual void request1DArray(unsigned int size);

	/**
	 * Return an uniformly distributed number on [0, 1).
	 * Undefined for purely deterministic samplers.
	 */
	virtual Float independent1D();

	/**
	 * Return an uniformly distributed 2D vector on [0, 1)x[0, 1).
	 * Undefined for purely deterministic samplers.
	 */
	virtual Point2 independent2D();

	/// Serialize this sampler to disk
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return a string description
	virtual std::string toString() const;

	/// Return the underlying random number generator
	inline Random *getRandom() { return m_random; }

	/**
	 * Update the current sample index, but without
	 * changing the RNG state. This is useful if the
	 * underlying random number generator has been used
	 * outside of this class
	 */
	inline void updateSampleIndex(uint64_t index) { m_sampleIndex = index; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~ReplayableSampler();
protected:
	ref<Random> m_initial, m_random;
};

MTS_NAMESPACE_END

#endif /* __REPLAYABLE_SAMPLER_H */
