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

#if !defined(__MEMPOOL_H)
#define __MEMPOOL_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * Basic memory pool -- allows repeated allocation & deallocation 
 * of objects of the same type, while attempting to keep them 
 * contiguous in memory and having only minimal interaction with the
 * underlying allocator.
 * \ingroup libcore
 */
template <typename T> class MemoryPool {
public:
	/// Create a new memory pool with aninitial set of 128 entries
	MemoryPool(size_t nEntries = 128) : m_size(0) {
		increaseCapacity(nEntries);
	}

	/// Destruct the memory pool and release all entries
	~MemoryPool() {
		for (size_t i=0; i<m_cleanup.size(); ++i)
			delete[] m_cleanup[i];
	}

	/// Acquire an entry
	inline T *alloc() {
		if (EXPECT_NOT_TAKEN(m_free.empty()))
			increaseCapacity();
		T *result = m_free.back();
		m_free.pop_back();
		return result;
	}

//	void ensureNotContained(T *ptr) {
//		if (std::find(m_free.begin(), m_free.end(), ptr) != m_free.end())
//			SLog(EError, "Memory pool inconsistency()");
//	}

	/// Release an entry
	inline void release(T *ptr) {
//		if (std::find(m_free.begin(), m_free.end(), ptr) != m_free.end())
//			SLog(EError, "Memory pool inconsistency in release()");
		m_free.push_back(ptr);
	}

	/// Return the total size of the memory pool
	inline size_t getSize() const {
		return m_size;
	}

	/// Check if every entry has been released
	bool isUnused() const {
		return m_free.size() == m_size;
	}
private:
	void increaseCapacity(size_t nEntries = 128) {
		T *ptr = new T[nEntries];
		for (size_t i=0; i<nEntries; ++i)
			m_free.push_back(&ptr[i]);
		m_cleanup.push_back(ptr);
		m_size += nEntries;
	}

private:
	std::vector<T *> m_free;
	std::vector<T *> m_cleanup;
	size_t m_size;
};

MTS_NAMESPACE_END

#endif /* __MEMPOOL_H */
