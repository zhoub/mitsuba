#include <mitsuba/core/fresolver.h>
#include <boost/algorithm/string.hpp>

#if defined(__WINDOWS__)
# include <windows.h>
# include <vector>
#endif

MTS_NAMESPACE_BEGIN

FileResolver::FileResolver() {
	m_paths.push_back(fs::current_path());
#if defined(__LINUX__)
	char exePathTemp[PATH_MAX];
	memset(exePathTemp, 0, PATH_MAX);
	if (readlink("/proc/self/exe", exePathTemp, PATH_MAX) != -1) {
		fs::path exePath(exePathTemp);

		/* Make sure that we're not running inside a Python interpreter */
		if (exePath.filename().string().find("python") == std::string::npos) {
			prependPath(exePath.parent_path());
			// Handle local installs: ~/local/bin/:~/local/share/mitsuba/*
			fs::path sharedDir = exePath.parent_path().parent_path()
				/ fs::path("share") / fs::path("mitsuba");
			if (fs::exists(sharedDir))
				prependPath(sharedDir);
		}
	} else {
		Log(EError, "Could not detect the executable path!");
	}
#elif defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	fs::path path = __mts_bundlepath();
	if (path.filename() != fs::path("Python.app"))
		prependPath(path);
	MTS_AUTORELEASE_END()
#elif defined(__WINDOWS__)
	std::vector<WCHAR> lpFilename(MAX_PATH);

	// Try to get the path with the default MAX_PATH length (260 chars)
	DWORD nSize = GetModuleFileNameW(NULL, &lpFilename[0], MAX_PATH);

	// Adjust the buffer size in case if was too short
	while (nSize == lpFilename.size()) {
		lpFilename.resize(nSize * 2);
		nSize = GetModuleFileNameW(NULL, &lpFilename[0], nSize);
	}

	// There is an error if and only if the function returns 0
	if (nSize != 0) {
		fs::path path(lpFilename);
		if (boost::to_lower_copy(path.filename().string()).find("python") == std::string::npos)
			prependPath(path.parent_path());
	} else {
		const std::string msg(lastErrorText());
		Log(EError, "Could not detect the executable path! (%s)", msg.c_str());
	}
#endif
}

FileResolver *FileResolver::clone() const {
	FileResolver *cloned = new FileResolver();
	cloned->m_paths = m_paths;
	return cloned;
}

void FileResolver::clear() {
	m_paths.clear();
}

void FileResolver::prependPath(const fs::path &path) {
	for (size_t i=0; i<m_paths.size(); ++i) {
		if (m_paths[i] == path)
			return;
	}
	m_paths.push_front(path);
}

void FileResolver::appendPath(const fs::path &path) {
	for (size_t i=0; i<m_paths.size(); ++i) {
		if (m_paths[i] == path)
			return;
	}
	m_paths.push_back(path);
}

fs::path FileResolver::resolve(const fs::path &path) const {
	for (size_t i=0; i<m_paths.size(); i++) {
		fs::path newPath = m_paths[i] / path;
		if (fs::exists(newPath))
			return newPath;
	}
	return path;
}

std::vector<fs::path> FileResolver::resolveAll(const fs::path &path) const {
	std::vector<fs::path> results;

	for (size_t i=0; i<m_paths.size(); i++) {
		fs::path newPath = m_paths[i] / path;
		if (fs::exists(newPath))
			results.push_back(newPath);
	}

	return results;
}

fs::path FileResolver::resolveAbsolute(const fs::path &path) const {
	return fs::absolute(resolve(path));
}

std::string FileResolver::toString() const {
	std::ostringstream oss;
	oss << "FileResolver[" << endl
		<< "  paths = {" << endl;
	for (size_t i=0; i<m_paths.size(); ++i) {
		oss << "    \"" << m_paths[i].string() << "\"";
		if (i+1 < m_paths.size())
			oss << ",";
		oss << endl;
	}
	oss << "  }" << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(FileResolver, false, Object)
MTS_NAMESPACE_END
