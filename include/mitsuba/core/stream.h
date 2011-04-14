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

#if !defined(__STREAM_H)
#define __STREAM_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract seekable stream class
 *
 * Specifies all functions to be implemented by stream 
 * subclasses and provides various convenience functions
 * based on them.
 *
 * All read<b>X</b>() and write<b>X</b>() methods support transparent
 * endianness conversion based on the endianness of the underlying 
 * system and the value passed to \ref setByteOrder().
 *
 * \sa FileStream, MemoryStream, SocketStream, 
 *     ConsoleStream, SSHStream, ZStream
 * \ingroup libcore
 */
class MTS_EXPORT_CORE Stream : public Object {
public:
	/// Defines the byte order to use in this Stream
	enum EByteOrder {
		EBigEndian = 0,                ///< PowerPC, SPARC, Motorola 68K
		ELittleEndian = 1,             ///< x86, x86_64
		ENetworkByteOrder = EBigEndian ///< Network byte order (an alias for big endian)
	};

	/**
	 * \brief Create a new stream. 
	 * 
	 * By default, it assumes the byte order of the 
	 * underlying system, i.e. no endianness conversion 
	 * is performed.
	 */
	Stream();
	
	/// Set the stream byte order
	void setByteOrder(EByteOrder byteOrder);

	/// Return the byte order of this stream
	inline EByteOrder getByteOrder() const { return m_byteOrder; }

	/// Return the byte order of the underlying machine
	inline EByteOrder getHostByteOrder() { return m_hostByteOrder; }

	/// Return a string representation
	virtual std::string toString() const;

	// ======================================================================
	/// @{ \name Abstract methods that need to be implemented by subclasses
	// ======================================================================

	/// Read a specified amount of data from the stream
	virtual void read(void *ptr, size_t size) = 0;

	/// Write a specified amount of data into the stream
	virtual void write(const void *ptr, size_t size) = 0;

	/// Seek to a position inside the stream
	virtual void setPos(size_t pos) = 0;

	/// Truncate the stream to a given size
	virtual void truncate(size_t size) = 0;

	/// Get the current position inside the stream
	virtual size_t getPos() const = 0;

	/// Return the size of the stream
	virtual size_t getSize() const = 0;

	/// Flush the stream's buffers
	virtual void flush() = 0;

	/// Can we write to the stream?
	virtual bool canWrite() const = 0;

	/// Can we read from the stream?
	virtual bool canRead() const = 0;

	/// @}
	// ======================================================================
	
	// ======================================================================
	/// @{ \name Convenience functions with automatic endianness conversion
	// ======================================================================

	/// Skip the given number of bytes
	void skip(size_t amount);

	/// Write a null-terminated string to the stream
	void writeString(const std::string &value);

	/// Write a string followed by a newline
	void writeLine(const std::string &value);

	/// Write a signed short (16 bit) to the stream
	void writeShort(short value);
	
	/// Write an array of signed shorts (16 bit) to the stream
	void writeShortArray(const short *values, size_t size);

	/// Write an unsigned short (16 bit) to the stream
	void writeUShort(unsigned short value);

	/// Write an array of unsigned shorts (16 bit) to the stream
	void writeUShortArray(const unsigned short *values, size_t size);

	/// Write a signed int (32 bit) to the stream
	void writeInt(int value);
	
	/// Write an array of signed ints (32 bit) to the stream
	void writeIntArray(const int *values, size_t size);

	/// Write an unsigned int (32 bit) to the stream
	void writeUInt(unsigned int value);
	
	/// Write an array of unsigned ints (32 bit) to the stream
	void writeUIntArray(const unsigned int *values, size_t size);

	/// Write a signed int (32 bit) to the stream
	void writeLong(int64_t value);

	/// Write an array of signed ints (64 bit) to the stream
	void writeLongArray(const int64_t *values, size_t size);

	/// Write an unsigned int (64 bit) to the stream
	void writeULong(uint64_t value);

	/// Write a size value to the stream
	void writeSize(size_t value) { writeULong((uint64_t) value); }

	/// Write an array of unsigned ints (64 bit) to the stream
	void writeULongArray(const uint64_t *values, size_t size);

	/// Write a signed character (8 bit) to the stream
	void writeChar(char value);

	/// Write an unsigned character (8 bit) to the stream
	void writeUChar(unsigned char value);

	/// Write a boolean (8 bit) to the stream
	inline void writeBool(bool value) { writeUChar(value); }
	
	/// Write a single-precision floating point number (32 bit) to the stream
	void writeSingle(float value);
	
	/// Write a single-precision floating point array (32 bit) to the stream
	void writeSingleArray(const float *data, size_t size);

	/// Write a double-precision floating point number (64 bit) to the stream
	void writeDouble(double value);
	
	/// Write a double-precision floating point array (64 bit) to the stream
	void writeDoubleArray(const double *data, size_t size);

	/// Write a floating point number (configured precision) to the stream
	inline void writeFloat(Float value) {
#ifdef SINGLE_PRECISION
		writeSingle(value);
#else
		writeDouble(value);
#endif
	}

	/// Write an array of floating point values (configured precision) to the stream
	inline void writeFloatArray(const Float *data, size_t size) {
#ifdef SINGLE_PRECISION
		writeSingleArray(data, size);
#else
		writeDoubleArray(data, size);
#endif
	}

	/// Return whether we are at the end of the stream
	bool isEOF() const;

	/// Read a line from the stream and return it as a string
	std::string readLine();

	/// Read a null-terminated string from the stream
	std::string readString();

	/// Read a signed short (16 bit) from the stream
	short readShort();
	
	/// Read an array of signed shorts (16 bit) from the stream
	void readShortArray(short *dest, size_t size);

	/// Read an unsigned short (16 bit) from the stream
	unsigned short readUShort();
	
	/// Read an array of unsigned shorts (16 bit) from the stream
	void readUShortArray(unsigned short *dest, size_t size);

	/// Read a signed int (32 bit) from the stream
	int readInt();
	
	/// Read an array of signed ints (32 bit) from the stream
	void readIntArray(int *dst, size_t size);

	/// Read an unsigned int (32 bit) from the stream
	unsigned int readUInt();
	
	/// Read an array of unsigned ints (32 bit) from the stream
	void readUIntArray(unsigned int *dest, size_t size);

	/// Read a signed int (64 bit) from the stream
	int64_t readLong();

	/// Read an array of signed ints (64 bit) from the stream
	void readLongArray(int64_t *dst, size_t size);

	/// Read an unsigned int (64 bit) from the stream
	uint64_t readULong();
	
	/// Read a size value from the stream
	size_t readSize() { return (size_t) readULong(); }

	/// Read an array of unsigned ints (64 bit) from the stream
	void readULongArray(uint64_t *dst, size_t size);

	/// Read a signed character (8 bit) from the stream
	char readChar();

	/// Read an unsigned character (8 bit) from the stream
	unsigned char readUChar();

	/// Read a boolean (8 bit) from the stream
	inline bool readBool() { return static_cast<bool> (readUChar()); }

	/// Read a single-precision floating point number (32 bit) from the stream
	float readSingle();

	/// Read a double-precision floating point array (64 bit) from the stream
	void readSingleArray(float *data, size_t size);

	/// Read a double-precision floating point number (64 bit) from the stream
	double readDouble();

	/// Read a double-precision floating point array (64 bit) from the stream
	void readDoubleArray(double *data, size_t size);

	/// Write a floating point number (configured precision) to the stream
	inline Float readFloat() {
#ifdef SINGLE_PRECISION
		return readSingle();
#else
		return readDouble();
#endif
	}

	/// Write an array of floating point values (configured precision) to the stream
	inline void readFloatArray(Float *data, size_t size) {
#ifdef SINGLE_PRECISION
		readSingleArray(data, size);
#else
		readDoubleArray(data, size);
#endif
	}

	/**
	 * \brief Copy content from this stream into another stream
	 * \param stream Destination stream
	 * \param numBytes 
	 * 		The number of bytes to copy. When -1 is specified,
	 * 		copying proceeds until the end of the source stream.
	 */
	void copyTo(Stream *stream, int64_t numBytes = -1);

	/**
	 * \brief Read an element from the stream (uses partial template 
	 * specialization to select a method appropriate to the data type)
	 */
	template <typename T> T readElement();

	/**
	 * \brief Write an element to the stream (uses partial template 
	 * specialization to select a method appropriate to the data type)
	 */
	template <typename T> void writeElement(T value);

	/**
	 * \brief Read an array from the stream (uses partial template 
	 * specialization to select a method appropriate to the data type)
	 */
	template <typename T> void readArray(T *array, size_t count);

	/**
	 * \brief Write an array to the stream (uses partial template 
	 * specialization to select a method appropriate to the data type)
	 */
	template <typename T> void writeArray(const T *array, size_t count);

	/// @}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Stream() { }
private:
	static EByteOrder m_hostByteOrder;
	EByteOrder m_byteOrder;
};


template <typename T> inline T Stream::readElement() {
	Log(EError, "Stream::readElement<T>: not implemented!");
}

template <typename T> inline void Stream::writeElement(T value) {
	Log(EError, "Stream::writeElement<T>: not implemented!");
}

template <typename T> inline void Stream::readArray(T *array, size_t count) {
	Log(EError, "Stream::readArray<T>: not implemented!");
}

template <typename T> inline void Stream::writeArray(const T *array, size_t count) {
	Log(EError, "Stream::writeArray<T>: not implemented!");
}

/// \cond
template <> inline float Stream::readElement() { return readSingle(); }
template <> inline double Stream::readElement() { return readDouble(); }
template <> inline char Stream::readElement() { return readChar(); }
template <> inline unsigned char Stream::readElement() { return readUChar(); }
template <> inline bool Stream::readElement() { return readBool(); }
template <> inline int16_t Stream::readElement() { return readShort(); }
template <> inline uint16_t Stream::readElement() { return readUShort(); }
template <> inline int32_t Stream::readElement() { return readInt(); }
template <> inline uint32_t Stream::readElement() { return readUInt(); }
template <> inline int64_t Stream::readElement() { return readLong(); }
template <> inline uint64_t Stream::readElement() { return readULong(); }
template <> inline std::string Stream::readElement() { return readString(); }

template <> inline void Stream::writeElement(float val) { return writeSingle(val); }
template <> inline void Stream::writeElement(double val) { return writeDouble(val); }
template <> inline void Stream::writeElement(char val) { return writeChar(val); }
template <> inline void Stream::writeElement(unsigned char val) { return writeUChar(val); }
template <> inline void Stream::writeElement(bool val) { return writeBool(val); }
template <> inline void Stream::writeElement(int16_t val) { return writeShort(val); }
template <> inline void Stream::writeElement(uint16_t val) { return writeUShort(val); }
template <> inline void Stream::writeElement(int32_t val) { return writeInt(val); }
template <> inline void Stream::writeElement(uint32_t val) { return writeUInt(val); }
template <> inline void Stream::writeElement(int64_t val) { return writeLong(val); }
template <> inline void Stream::writeElement(uint64_t val) { return writeULong(val); }
template <> inline void Stream::writeElement(const std::string &val) { return writeString(val); }

template <> inline void Stream::readArray(float *array, size_t count) { return readSingleArray(array, count); }
template <> inline void Stream::readArray(double *array, size_t count) { return readDoubleArray(array, count); }
template <> inline void Stream::readArray(char *array, size_t count) { return read((uint8_t *) array, count); }
template <> inline void Stream::readArray(unsigned char *array, size_t count) { return read((uint8_t *) array, count); }
template <> inline void Stream::readArray(bool *array, size_t count) { return read((uint8_t *) array, count); }
template <> inline void Stream::readArray(int16_t *array, size_t count) { return readShortArray(array, count); }
template <> inline void Stream::readArray(uint16_t *array, size_t count) { return readUShortArray(array, count); }
template <> inline void Stream::readArray(int32_t *array, size_t count) { return readIntArray(array, count); }
template <> inline void Stream::readArray(uint32_t *array, size_t count) { return readUIntArray(array, count); }
template <> inline void Stream::readArray(int64_t *array, size_t count) { return readLongArray(array, count); }
template <> inline void Stream::readArray(uint64_t *array, size_t count) { return readULongArray(array, count); }

template <> inline void Stream::writeArray(const float *array, size_t count) { return writeSingleArray(array, count); }
template <> inline void Stream::writeArray(const double *array, size_t count) { return writeDoubleArray(array, count); }
template <> inline void Stream::writeArray(const char *array, size_t count) { return write((uint8_t *) array, count); }
template <> inline void Stream::writeArray(const unsigned char *array, size_t count) { return write((uint8_t *) array, count); }
template <> inline void Stream::writeArray(const bool *array, size_t count) { return write((uint8_t *) array, count); }
template <> inline void Stream::writeArray(const int16_t *array, size_t count) { return writeShortArray(array, count); }
template <> inline void Stream::writeArray(const uint16_t *array, size_t count) { return writeUShortArray(array, count); }
template <> inline void Stream::writeArray(const int32_t *array, size_t count) { return writeIntArray(array, count); }
template <> inline void Stream::writeArray(const uint32_t *array, size_t count) { return writeUIntArray(array, count); }
template <> inline void Stream::writeArray(const int64_t *array, size_t count) { return writeLongArray(array, count); }
template <> inline void Stream::writeArray(const uint64_t *array, size_t count) { return writeULongArray(array, count); }


/// \endcond

MTS_NAMESPACE_END

#endif /* __STREAM_H */
