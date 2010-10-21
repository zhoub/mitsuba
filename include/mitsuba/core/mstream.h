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

#if !defined(__MSTREAM_H)
#define __MSTREAM_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple memory buffer-based stream with automatic memory management 
 * 
 * The underlying memory storage of this implementation dynamically expands 
 * as data is written to the stream.
 */
class MTS_EXPORT_CORE MemoryStream : public Stream {
public:
	/// Create a new memory stream
	MemoryStream(size_t initialSize = 512);

	/// Create a memory stream, which operates on a pre-allocated buffer
	MemoryStream(uint8_t *ptr, size_t size);

	/// Return the underlying data
	inline uint8_t *getData() { return m_data; }
	
	/// Return the underlying data
	inline const uint8_t *getData() const { return m_data; }

	/// Set size and position to zero without changing the underlying buffer
	void reset();

	/* Stream implementation */
	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void setPos(size_t pos);
	size_t getPos() const;
	size_t getSize() const;
	void truncate(size_t size);
	void flush();
	bool canWrite() const;
	bool canRead() const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	void resize(size_t newSize);

	// \brief Virtual destructor
	virtual ~MemoryStream();
protected:
	size_t m_capacity;
	size_t m_size;
	size_t m_pos;
	bool m_ownsBuffer;
	uint8_t *m_data;
};

MTS_NAMESPACE_END

#endif /* __MSTREAM_H */
