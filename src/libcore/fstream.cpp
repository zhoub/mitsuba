/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/core/fstream.h>
#include <cerrno>

#if !defined(WIN32)
# include <unistd.h>
#else
# include <windows.h>
#endif

MTS_NAMESPACE_BEGIN

struct FileStream::FileStreamPrivate
{
#if defined(WIN32)
	HANDLE file;
#else
	FILE* file;
#endif
	bool write;
	bool read;
	FileStream::EFileMode mode;
	fs::path path;

	FileStreamPrivate() : file(NULL) {}
};

FileStream::FileStream()
 : d(new FileStreamPrivate) {
}

FileStream::FileStream(const fs::path &path, EFileMode mode)
 : d(new FileStreamPrivate) {
	open(path, mode);
}


FileStream::~FileStream() {
	if (d->file != 0)
		close();
}

const fs::path& FileStream::getPath() const {
	return d->path;
}

std::string FileStream::toString() const {
	std::ostringstream oss;
	oss << "FileStream[" << Stream::toString()
		<< ", path=\"" << d->path.string()
		<< "\", mode=" << d->mode << "]";
	return oss.str();
}

void FileStream::open(const fs::path &path, EFileMode mode) {
	AssertEx(d->file == 0, "A file has already been opened using this stream");

	Log(ETrace, "Opening \"%s\"", path.string().c_str());

	d->path = path;
	d->mode = mode;
	d->write = true;
	d->read = true;

#ifdef WIN32
	DWORD dwDesiredAccess = GENERIC_READ;
	DWORD dwCreationDisposition = OPEN_EXISTING;

	switch (d->mode) {
	case EReadOnly:
		d->write = false;
		break;
	case EReadWrite:
		dwDesiredAccess |= GENERIC_WRITE;
		break;
	case ETruncWrite:
		d->read = false;
		dwDesiredAccess = GENERIC_WRITE;
		dwCreationDisposition = CREATE_ALWAYS;
		break;
	case ETruncReadWrite:
		dwDesiredAccess |= GENERIC_WRITE;
		dwCreationDisposition = CREATE_ALWAYS;
		break;
	case EAppendWrite:
		d->read = false;
		dwDesiredAccess = GENERIC_WRITE;
		break;
	case EAppendReadWrite:
		dwDesiredAccess |= GENERIC_WRITE;
		break;
	default: 
		Log(EError, "Unknown file mode");
		break;
	}

	d->file = CreateFile(path.string().c_str(), dwDesiredAccess, 
		FILE_SHARE_WRITE | FILE_SHARE_READ, 0, 
		dwCreationDisposition, FILE_ATTRIBUTE_NORMAL, 0);

	if (d->file == INVALID_HANDLE_VALUE)
		Log(EError, "Error while trying to open file \"%s\": %s", 
			d->path.string().c_str(), lastErrorText().c_str());
	
	if (d->mode == EAppendWrite || d->mode == EAppendReadWrite)
		seek(getSize());
#else
	const char *modeString = NULL;

	switch (d->mode) {
	case EReadOnly:
		modeString = "rb";
		d->write = false;
		break;
	case EReadWrite:
		modeString = "rb+";
		break;
	case ETruncWrite:
		modeString = "wb";
		d->read = false;
		break;
	case ETruncReadWrite:
		modeString = "wb+";
		break;
	case EAppendWrite:
		modeString = "ab";
		d->read = false;
		break;
	case EAppendReadWrite:
		modeString = "ab+";
		break;
	default:
		Log(EError, "Unknown file mode");
		break;
	};

	d->file = fopen(d->path.string().c_str(), modeString);

	if (d->file == NULL) {
		Log(EError, "Error while trying to open file \"%s\": %s", 
			d->path.string().c_str(), strerror(errno));
	}
#endif
}

void FileStream::close() {
	AssertEx(d->file != 0, "No file is currently open");
	Log(ETrace, "Closing \"%s\"", d->path.string().c_str());

#ifdef WIN32
	if (!CloseHandle(d->file)) {
		Log(EError, "Error while trying to close file \"%s\": %s", 
			d->path.string().c_str(), lastErrorText().c_str());
	}
#else
	if (fclose(d->file)) {
		Log(EError, "Error while trying to close file \"%s\": %s", 
			d->path.string().c_str(), strerror(errno));
	}
#endif
	d->file = 0;
}


void FileStream::remove() {
	close();
	Log(EDebug, "Removing \"%s\"", 	d->path.string().c_str());

	fs::remove(d->path);
}

void FileStream::seek(size_t pos) {
	AssertEx(d->file != 0, "No file is currently open");
	
#ifdef WIN32
	LARGE_INTEGER fpos;
	fpos.QuadPart = pos;
	if (SetFilePointerEx(d->file, fpos, 0, FILE_BEGIN) == INVALID_SET_FILE_POINTER) {
		Log(EError, "Error while trying to seek to position %i in file \"%s\": %s", 
			pos, d->path.string().c_str(), lastErrorText().c_str());
	}
#else
	if (fseek(d->file, pos, SEEK_SET)) {
		Log(EError, "Error while trying to seek to position %i in file \"%s\": %s", 
			pos, d->path.string().c_str(), strerror(errno));
	}
#endif
}

size_t FileStream::getPos() const {
	AssertEx(d->file != 0, "No file is currently open");
#ifdef WIN32
	DWORD pos = SetFilePointer(d->file, 0, 0, FILE_CURRENT);
	if (pos == INVALID_SET_FILE_POINTER) {
		Log(EError, "Error while looking up the position in file \"%s\": %s",
			d->path.string().c_str(), lastErrorText().c_str());
	}
	return (size_t) pos;
#else
	long pos;
	pos = ftell(d->file);
	if (pos == -1) {
		Log(EError, "Error while looking up the position in file \"%s\": %s", 
			d->path.string().c_str(), strerror(errno));
	}
	return (size_t) pos;
#endif
}

size_t FileStream::getSize() const {
	AssertEx(d->file != 0, "No file is currently open");

#ifdef WIN32
	LARGE_INTEGER result;
	if (GetFileSizeEx(d->file, &result) == 0) {
		Log(EError, "Error while getting the file size of \"%s\": %s",
			d->path.string().c_str(), lastErrorText().c_str());
	}
	return (size_t) result.QuadPart;
#else
	size_t size, tmp;
	
	tmp = getPos();
	if (fseek(d->file, 0, SEEK_END)) {
		Log(EError, "Error while seeking within \"%s\": %s",
			d->path.string().c_str(), strerror(errno));
	}
	size = getPos();
	if (fseek(d->file, tmp, SEEK_SET)) {
		Log(EError, "Error while seeking within \"%s\": %s",
			d->path.string().c_str(), strerror(errno));
	}
	return size;
#endif
}

void FileStream::truncate(size_t size) {
	AssertEx(d->file != 0, "No file is currently open");
	AssertEx(d->write, "File is not open with write access");

	size_t pos = getPos();
	if (pos > size) 
		pos = size;

#ifdef WIN32
	seek(size);
	if (!SetEndOfFile(d->file)) {
		Log(EError, "Error while truncating file \"%s\": %s",
			d->path.string().c_str(), lastErrorText().c_str());
	}
#else
	seek(pos);
	flush();

	if (ftruncate(fileno(d->file), size)) {
		Log(EError, "Error while truncating file \"%s\": %s",
			d->path.string().c_str(), strerror(errno));
	}
#endif
	seek(pos);
}

void FileStream::flush() {
	AssertEx(d->file != 0, "No file is currently open");
	AssertEx(d->write, "File is not open with write access");
#ifdef WIN32
	if (!FlushFileBuffers(d->file)) {
		Log(EError, "Error while flusing the buffers of \"%s\": %s",
			d->path.string().c_str(), lastErrorText().c_str());
	}
#else
	if (fflush(d->file) != 0) {
		Log(EError, "Error while flusing the buffers of \"%s\": %s",
			d->path.string().c_str(), strerror(errno));
	}
#endif
}

void FileStream::read(void *pPtr, size_t size) {
	AssertEx(d->file != 0, "No file is currently open");
	AssertEx(d->read, "File is not open with read access");
	
	if (size == 0)
		return;
#ifdef WIN32
	DWORD lpNumberOfBytesRead;
	if (!ReadFile(d->file, pPtr, (DWORD) size, &lpNumberOfBytesRead, 0)) {
		Log(EError, "Error while reading from file \"%s\": %s",
			d->path.string().c_str(), lastErrorText().c_str());
	}
	if (lpNumberOfBytesRead != (DWORD) size) 
		throw EOFException(formatString("Read less data than expected (%i bytes required) "
			"from file \"%s\"", size, d->path.string().c_str()), (size_t) lpNumberOfBytesRead);
#else
	size_t bytesRead;
	if ((bytesRead = fread(pPtr, 1, size, d->file)) != size) {
		if (ferror(d->file) != 0) {
			Log(EError, "Error while reading from file \"%s\": %s",
				d->path.string().c_str(), strerror(errno));
		}
		throw EOFException(formatString("Read less data than expected (%i bytes required) "
			"from file \"%s\"", size, d->path.string().c_str()), bytesRead);
	}
#endif
}

void FileStream::write(const void *pPtr, size_t size) {
	AssertEx(d->file != 0, "No file is currently open");
	AssertEx(d->write, "File is not open with write access");

	if (size == 0)
		return;

#ifdef WIN32
	DWORD lpNumberOfBytesWritten;
	if (!WriteFile(d->file, pPtr, (DWORD) size, &lpNumberOfBytesWritten, 0)) {
		Log(EError, "Error while writing to file \"%s\": %s",
			d->path.string().c_str(), lastErrorText().c_str());
	}
	if (lpNumberOfBytesWritten != (DWORD) size) 
		throw EOFException(formatString("Wrote less data than expected (%i bytes required) "
			"to file \"%s\"", size, d->path.string().c_str()), (size_t) lpNumberOfBytesWritten);
#else
	size_t bytesWritten;
	if ((bytesWritten = fwrite(pPtr, 1, size, d->file)) != size) {
		if (ferror(d->file))
			Log(EError, "Error while writing to file \"%s\": %s",
				d->path.string().c_str(), strerror(errno));
		throw EOFException(formatString("Wrote less data than expected (%i bytes required) "
			"to file \"%s\"", size, d->path.string().c_str()), bytesWritten);
	}
#endif
}

bool FileStream::canRead() const {
	AssertEx(d->file != 0, "No file is currently open");
	return d->read;
}

bool FileStream::canWrite() const {
	AssertEx(d->file != 0, "No file is currently open");
	return d->write;
}

MTS_IMPLEMENT_CLASS(FileStream, false, Stream)
MTS_NAMESPACE_END
