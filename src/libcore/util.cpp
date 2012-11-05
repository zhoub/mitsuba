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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/util.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/sse.h>
#include <boost/bind.hpp>
#include <stdarg.h>
#include <iomanip>
#include <errno.h>

#if defined(__OSX__)
#include <sys/sysctl.h>
#elif defined(WIN32)
#include <direct.h>
#else
#include <malloc.h>
#endif

#if defined(WIN32)
# include <windows.h>
# include <winsock2.h>
# include <ws2tcpip.h>
#else
# include <sys/types.h>
# include <sys/socket.h>
# include <netdb.h>
# include <fenv.h>
#endif

// SSE is not enabled in general when using double precision, however it is
// required in OS X for FP exception handling
#if defined(__OSX__) && !defined(MTS_SSE)
#include <xmmintrin.h>
#undef enable_fpexcept_sse
#undef query_fpexcept_sse
#undef disable_fpexcept_sse

namespace {
inline void enable_fpexcept_sse() {
	_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
		~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO));
}
inline unsigned int query_fpexcept_sse() {
	return (~_MM_GET_EXCEPTION_MASK() &
		(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO));
}
inline void disable_fpexcept_sse() {
	_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() |
		_MM_MASK_INVALID | _MM_MASK_DIV_ZERO);
}
} // namespace

#endif

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  General utility functions
// -----------------------------------------------------------------------

std::vector<std::string> tokenize(const std::string &string, const std::string &delim) {
	std::string::size_type lastPos = string.find_first_not_of(delim, 0);
	std::string::size_type pos = string.find_first_of(delim, lastPos);
	std::vector<std::string> tokens;

	while (std::string::npos != pos || std::string::npos != lastPos) {
		tokens.push_back(string.substr(lastPos, pos - lastPos));
		lastPos = string.find_first_not_of(delim, pos);
		pos = string.find_first_of(delim, lastPos);
	}

	return tokens;
}

std::string trim(const std::string& str) {
	std::string::size_type
		start = str.find_first_not_of(" \t\r\n"),
		end = str.find_last_not_of(" \t\r\n");

	return str.substr(start == std::string::npos ? 0 : start,
			end == std::string::npos ? str.length() - 1 : end - start + 1);
}

std::string indent(const std::string &string, int amount) {
	/* This could probably be done faster (is not
	   really speed-critical though) */
	std::istringstream iss(string);
	std::ostringstream oss;
	std::string str;
	bool firstLine = true;
	while (!iss.eof()) {
		std::getline(iss, str);
		if (!firstLine) {
			for (int i=0; i<amount; ++i)
				oss << "  ";
		}
		oss << str;
		if (!iss.eof())
			oss << endl;
		firstLine = false;
	}
	return oss.str();
}

std::string memString(size_t size) {
	Float value = (Float) size;
	const char *prefixes[] = {
		"B", "KiB", "MiB", "GiB", "TiB", "PiB"
	};
	int prefix = 0;
	while (prefix < 5 && value > 1024.0f) {
		value /= 1024.0f; ++prefix;
	}
	return formatString(prefix == 0 ?
			"%.0f %s" : "%.2f %s", value, prefixes[prefix]);
}

void * __restrict allocAligned(size_t size) {
#if defined(WIN32)
	return _aligned_malloc(size, L1_CACHE_LINE_SIZE);
#elif defined(__OSX__)
	/* OSX malloc already returns 16-byte aligned data suitable
	   for AltiVec and SSE computations */
	return malloc(size);
#else
	return memalign(L1_CACHE_LINE_SIZE, size);
#endif
}

void freeAligned(void *ptr) {
#if defined(WIN32)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

int getCoreCount() {
#if defined(WIN32)
	SYSTEM_INFO sys_info;
	GetSystemInfo(&sys_info);
	return sys_info.dwNumberOfProcessors;
#elif defined(__OSX__)
	int nprocs;
	size_t nprocsSize = sizeof(int);
	if (sysctlbyname("hw.activecpu", &nprocs, &nprocsSize, NULL, 0))
		SLog(EError, "Could not detect the number of processors!");
	return (int) nprocs;
#else
	return sysconf(_SC_NPROCESSORS_CONF);
#endif
}

#if defined(WIN32)
std::string lastErrorText() {
	DWORD errCode = GetLastError();
	char *errorText = NULL;
	if (!FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER
		| FORMAT_MESSAGE_FROM_SYSTEM
		| FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		errCode,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &errorText,
		0,
		NULL)) {
		return "Internal error while looking up an error code";
	}
	std::string result(errorText);
	LocalFree(errorText);
	return result;
}
#endif

bool enableFPExceptions() {
	bool exceptionsWereEnabled = false;
#if defined(WIN32)
	_clearfp();
	uint32_t cw = _controlfp(0, 0);
	exceptionsWereEnabled = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	cw &= ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	_controlfp(cw, _MCW_EM);
#elif defined(__OSX__)
	exceptionsWereEnabled = query_fpexcept_sse() != 0;
#else
	exceptionsWereEnabled =
		fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
	enable_fpexcept_sse();
	return exceptionsWereEnabled;
}

bool disableFPExceptions() {
	bool exceptionsWereEnabled = false;
#if defined(WIN32)
	_clearfp();
	uint32_t cw = _controlfp(0, 0);
	exceptionsWereEnabled = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	cw |= _EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW;
	_controlfp(cw, _MCW_EM);
#elif defined(__OSX__)
	exceptionsWereEnabled = query_fpexcept_sse() != 0;
#else
	exceptionsWereEnabled =
		fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
	fedisableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
	disable_fpexcept_sse();
	return exceptionsWereEnabled;
}

void restoreFPExceptions(bool oldState) {
	bool currentState;
#if defined(WIN32)
	uint32_t cw = _controlfp(0, 0);
	currentState = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
#elif defined(__OSX__)
	currentState = query_fpexcept_sse() != 0;
#else
	currentState = fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
	if (oldState != currentState) {
		if (oldState)
			enableFPExceptions();
		else
			disableFPExceptions();
	}
}

std::string getHostName() {
	char hostName[128];
	if (gethostname(hostName, sizeof(hostName)) != 0)
#if defined(WIN32)
		SLog(EError, "Could not retrieve the computer's host name: %s!",
			lastErrorText().c_str());
#else
		SLog(EError, "Could not retrieve the computer's host name : %s!",
			strerror(errno));
#endif
	return hostName;
}

std::string getFQDN() {
	struct addrinfo *addrInfo = NULL, hints;
	memset(&hints, 0, sizeof(addrinfo));
	// Only look for IPv4 addresses -- perhaps revisit this later
	hints.ai_family = AF_INET; // AF_UNSPEC
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;

	int retVal = getaddrinfo(getHostName().c_str(), NULL, &hints, &addrInfo);
	if (addrInfo == NULL || retVal != 0) {
		SLog(EWarn, "Could not retrieve the computer's fully "
			"qualified domain name: could not resolve host address \"%s\"!",
			getHostName().c_str());
		return getHostName();
	}

	char fqdn[NI_MAXHOST];
	retVal = getnameinfo(addrInfo->ai_addr, sizeof(struct sockaddr_in),
		fqdn, NI_MAXHOST, NULL, 0, 0);
	if (retVal != 0) {
		freeaddrinfo(addrInfo);
#if defined(WIN32)
		SLog(EWarn, "Could not retrieve the computer's fully "
			"qualified domain name: error %i!", WSAGetLastError());
#else
		SLog(EWarn, "Could not retrieve the computer's fully "
			"qualified domain name: error %i!", gai_strerror(retVal));
#endif
		return getHostName();
	}

	freeaddrinfo(addrInfo);

	return fqdn;
}

Float log2(Float value) {
	const Float invLn2 = (Float) 1.0f / math::fastlog((Float) 2.0f);
	return math::fastlog(value) * invLn2;
}

std::string formatString(const char *fmt, ...) {
	char tmp[512];
	va_list iterator;

#if defined(WIN32)
	va_start(iterator, fmt);
	size_t size = _vscprintf(fmt, iterator) + 1;

	if (size >= sizeof(tmp)) {
		char *dest = new char[size];
		vsnprintf_s(dest, size, size-1, fmt, iterator);
		va_end(iterator);
		std::string result(dest);
		delete[] dest;
		return result;
	}

	vsnprintf_s(tmp, size, size-1, fmt, iterator);
	va_end(iterator);
#else
	va_start(iterator, fmt);
	size_t size = vsnprintf(tmp, sizeof(tmp), fmt, iterator);
	va_end(iterator);

	if (size >= sizeof(tmp)) {
		/* Overflow! -- dynamically allocate memory */
		char *dest = new char[size+1];
		va_start(iterator, fmt);
		vsnprintf(dest, size+1, fmt, iterator);
		va_end(iterator);

		std::string result(dest);
		delete[] dest;
		return result;
	}
#endif

	return std::string(tmp);
}

int log2i(uint32_t value) {
	int r = 0;
	while ((value >> r) != 0)
		r++;
	return r-1;
}

int log2i(uint64_t value) {
	int r = 0;
	while ((value >> r) != 0)
		r++;
	return r-1;
}

/* Fast rounding & power-of-two test algorithms from PBRT */
uint32_t roundToPowerOfTwo(uint32_t i) {
	i--;
	i |= i >> 1; i |= i >> 2;
	i |= i >> 4; i |= i >> 8;
	i |= i >> 16;
	return i+1;
}

uint64_t roundToPowerOfTwo(uint64_t i) {
	i--;
	i |= i >> 1;  i |= i >> 2;
	i |= i >> 4;  i |= i >> 8;
	i |= i >> 16; i |= i >> 32;
	return i+1;
}

// -----------------------------------------------------------------------
//  Numerical utility functions
// -----------------------------------------------------------------------

bool solveQuadratic(Float a, Float b, Float c, Float &x0, Float &x1) {
	/* Linear case */
	if (a == 0) {
		if (b != 0) {
			x0 = x1 = -c / b;
			return true;
		}
		return false;
	}

	Float discrim = b*b - 4.0f*a*c;

	/* Leave if there is no solution */
	if (discrim < 0)
		return false;

	Float temp, sqrtDiscrim = std::sqrt(discrim);

	/* Numerically stable version of (-b (+/-) sqrtDiscrim) / (2 * a)
	 *
	 * Based on the observation that one solution is always
	 * accurate while the other is not. Finds the solution of
	 * greater magnitude which does not suffer from loss of
	 * precision and then uses the identity x1 * x2 = c / a
	 */
	if (b < 0)
		temp = -0.5f * (b - sqrtDiscrim);
	else
		temp = -0.5f * (b + sqrtDiscrim);

	x0 = temp / a;
	x1 = c / temp;

	/* Return the results so that x0 < x1 */
	if (x0 > x1)
		std::swap(x0, x1);

	return true;
}

bool solveQuadraticDouble(double a, double b, double c, double &x0, double &x1) {
	/* Linear case */
	if (a == 0) {
		if (b != 0) {
			x0 = x1 = -c / b;
			return true;
		}
		return false;
	}

	double discrim = b*b - 4.0f*a*c;

	/* Leave if there is no solution */
	if (discrim < 0)
		return false;

	double temp, sqrtDiscrim = std::sqrt(discrim);

	/* Numerically stable version of (-b (+/-) sqrtDiscrim) / (2 * a)
	 *
	 * Based on the observation that one solution is always
	 * accurate while the other is not. Finds the solution of
	 * greater magnitude which does not suffer from loss of
	 * precision and then uses the identity x1 * x2 = c / a
	 */
	if (b < 0)
		temp = -0.5f * (b - sqrtDiscrim);
	else
		temp = -0.5f * (b + sqrtDiscrim);

	x0 = temp / a;
	x1 = c / temp;

	/* Return the results so that x0 < x1 */
	if (x0 > x1)
		std::swap(x0, x1);

	return true;
}

bool solveLinearSystem2x2(const Float a[2][2], const Float b[2], Float x[2]) {
	Float det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

	if (det == 0)
		return false;

	Float inverse = (Float) 1.0f / det;

	x[0] = (a[1][1] * b[0] - a[0][1] * b[1]) * inverse;
	x[1] = (a[0][0] * b[1] - a[1][0] * b[0]) * inverse;

	return true;
}

Float interpCubic1D(Float x, const Float *data, Float min, Float max, size_t size, bool extrapolate) {
	/* Give up when given an out-of-range or NaN argument */
	if (!(x >= min && x <= max) && !extrapolate)
		return 0.0f;

	/* Transform 'x' so that knots lie at integer positions */
	Float t = ((x - min) * (size - 1)) / (max - min);

	/* Find the index of the left knot in the queried subinterval, be
	   robust to cases where 't' lies exactly on the right endpoint */
	size_t k = std::max((size_t) 0, std::min((size_t) t, size - 2));

	Float f0  = data[k],
		  f1  = data[k+1],
		  d0, d1;

	/* Approximate the derivatives */
	if (k > 0)
		d0 = 0.5f * (data[k+1] - data[k-1]);
	else
		d0 = data[k+1] - data[k];

	if (k + 2 < size)
		d1 = 0.5f * (data[k+2] - data[k]);
	else
		d1 = data[k+1] - data[k];

	/* Compute the relative position within the interval */
	t = t - (Float) k;

	Float t2 = t*t, t3 = t2*t;

	return
		( 2*t3 - 3*t2 + 1) * f0 +
		(-2*t3 + 3*t2)     * f1 +
		(   t3 - 2*t2 + t) * d0 +
		(   t3 - t2)       * d1;
}

Float interpCubic1DIrregular(Float x, const Float *nodes, const Float *data, size_t size, bool extrapolate) {
	/* Give up when given an out-of-range or NaN argument */
	if (!(x >= nodes[0] && x <= nodes[size-1]) && !extrapolate)
		return 0.0f;

	size_t k = (size_t) std::max((ptrdiff_t) 0, std::min((ptrdiff_t) size - 2,
				std::lower_bound(nodes, nodes + size, x) - nodes - 1));

	Float f0       = data[k],
		  f1       = data[k+1],
		  width    = nodes[k+1] - nodes[k],
		  invWidth = 1.0f / width,
		  d0, d1;

	/* Approximate the derivatives */
	if (k > 0)
		d0 = (f1 - data[k-1]) / (nodes[k+1] - nodes[k-1]);
	else
		d0 = (f1 - f0) * invWidth;

	if (k + 2 < size)
		d1 = (data[k+2] - f0) / (nodes[k+2] - nodes[k]);
	else
		d1 = (f1 - f0) * invWidth;

	Float t = (x - nodes[k]) * invWidth;
	Float t2 = t*t, t3 = t2*t;

	return
	    ( 2*t3 - 3*t2 + 1) * f0 +
	    (-2*t3 + 3*t2)     * f1 +
	   ((   t3 - 2*t2 + t) * d0 +
	    (   t3 - t2)       * d1) * width;
}


Float interpCubic2D(const Point2 &p, const Float *data,
		const Point2 &min, const Point2 &max, const Size2 &size, bool extrapolate) {
	Float knotWeights[2][4];
	Size2 knot;

	/* Compute interpolation weights separately for each dimension */
	for (int dim=0; dim<2; ++dim) {
		Float *weights = knotWeights[dim];
		/* Give up when given an out-of-range or NaN argument */
		if (!(p[dim] >= min[dim] && p[dim] <= max[dim]) && !extrapolate)
			return 0.0f;

		/* Transform 'p' so that knots lie at integer positions */
		Float t = ((p[dim] - min[dim]) * (size[dim] - 1))
			/ (max[dim]-min[dim]);

		/* Find the index of the left knot in the queried subinterval, be
		   robust to cases where 't' lies exactly on the right endpoint */
		knot[dim] = std::min((size_t) t, size[dim] - 2);

		/* Compute the relative position within the interval */
		t = t - (Float) knot[dim];

		/* Compute node weights */
		Float t2 = t*t, t3 = t2*t;
		weights[0] = 0.0f;
		weights[1] = 2*t3 - 3*t2 + 1;
		weights[2] = -2*t3 + 3*t2;
		weights[3] = 0.0f;

		/* Derivative weights */
		Float d0 = t3 - 2*t2 + t,
			  d1 = t3 - t2;

		/* Turn derivative weights into node weights using
		   an appropriate chosen finite differences stencil */
		if (knot[dim] > 0) {
			weights[2] +=  0.5f * d0;
			weights[0] -=  0.5f * d0;
		} else {
			weights[2] += d0;
			weights[1] -= d0;
		}

		if (knot[dim] + 2 < size[dim]) {
			weights[3] += 0.5f * d1;
			weights[1] -= 0.5f * d1;
		} else {
			weights[2] += d1;
			weights[1] -= d1;
		}
	}

	Float result = 0.0f;
	for (int y=-1; y<=2; ++y) {
		Float wy = knotWeights[1][y+1];
		for (int x=-1; x<=2; ++x) {
			Float wxy = knotWeights[0][x+1] * wy;

			if (wxy == 0)
				continue;

			size_t pos = (knot[1] + y) * size[0] + knot[0] + x;

			result += data[pos] * wxy;
		}
	}
	return result;
}

Float interpCubic2DIrregular(const Point2 &p, const Float **nodes_,
			const Float *data, const Size2 &size, bool extrapolate) {
	Float knotWeights[2][4];
	Size2 knot;

	/* Compute interpolation weights separately for each dimension */
	for (int dim=0; dim<2; ++dim) {
		const Float *nodes = nodes_[dim];
		Float *weights = knotWeights[dim];

		/* Give up when given an out-of-range or NaN argument */
		if (!(p[dim] >= nodes[0] && p[dim] <= nodes[size[dim]-1]) && !extrapolate)
			return 0.0f;

		/* Find the index of the left knot in the queried subinterval, be
		   robust to cases where 't' lies exactly on the right endpoint */
		size_t k = (size_t) std::max((ptrdiff_t) 0, std::min((ptrdiff_t) size[dim] - 2,
			std::lower_bound(nodes, nodes + size[dim], p[dim]) - nodes - 1));
		knot[dim] = k;

		Float width = nodes[k+1] - nodes[k], invWidth = 1 / width;

		/* Compute the relative position within the interval */
		Float t = (p[dim] - nodes[k]) * invWidth,
			  t2 = t*t, t3 = t2*t;

		/* Compute node weights */
		weights[0] = 0.0f;
		weights[1] = 2*t3 - 3*t2 + 1;
		weights[2] = -2*t3 + 3*t2;
		weights[3] = 0.0f;

		/* Derivative weights */
		Float d0 = (t3 - 2*t2 + t) * width,
			  d1 = (t3 - t2) * width;

		/* Turn derivative weights into node weights using
		   an appropriate chosen finite differences stencil */
		if (k > 0) {
			Float factor = 1 / (nodes[k+1]-nodes[k-1]);
			weights[2] += d0 * factor;
			weights[0] -= d0 * factor;
		} else {
			weights[2] += d0 * invWidth;
			weights[1] -= d0 * invWidth;
		}

		if (k + 2 < size[dim]) {
			Float factor = 1 / (nodes[k+2]-nodes[k]);
			weights[3] += d1 * factor;
			weights[1] -= d1 * factor;
		} else {
			weights[2] += d1 * invWidth;
			weights[1] -= d1 * invWidth;
		}
	}

	Float result = 0.0f;
	for (int y=-1; y<=2; ++y) {
		Float wy = knotWeights[1][y+1];
		for (int x=-1; x<=2; ++x) {
			Float wxy = knotWeights[0][x+1] * wy;

			if (wxy == 0)
				continue;

			size_t pos = (knot[1] + y) * size[0] + knot[0] + x;

			result += data[pos] * wxy;
		}
	}
	return result;
}

Float interpCubic3D(const Point3 &p, const Float *data,
		const Point3 &min, const Point3 &max, const Size3 &size, bool extrapolate) {
	Float knotWeights[3][4];
	Size3 knot;

	/* Compute interpolation weights separately for each dimension */
	for (int dim=0; dim<3; ++dim) {
		Float *weights = knotWeights[dim];
		/* Give up when given an out-of-range or NaN argument */
		if (!(p[dim] >= min[dim] && p[dim] <= max[dim]) && !extrapolate)
			return 0.0f;

		/* Transform 'p' so that knots lie at integer positions */
		Float t = ((p[dim] - min[dim]) * (size[dim] - 1))
			/ (max[dim]-min[dim]);

		/* Find the index of the left knot in the queried subinterval, be
		   robust to cases where 't' lies exactly on the right endpoint */
		knot[dim] = std::min((size_t) t, size[dim] - 2);

		/* Compute the relative position within the interval */
		t = t - (Float) knot[dim];

		/* Compute node weights */
		Float t2 = t*t, t3 = t2*t;
		weights[0] = 0.0f;
		weights[1] = 2*t3 - 3*t2 + 1;
		weights[2] = -2*t3 + 3*t2;
		weights[3] = 0.0f;

		/* Derivative weights */
		Float d0 = t3 - 2*t2 + t,
			  d1 = t3 - t2;

		/* Turn derivative weights into node weights using
		   an appropriate chosen finite differences stencil */
		if (knot[dim] > 0) {
			weights[2] +=  0.5f * d0;
			weights[0] -=  0.5f * d0;
		} else {
			weights[2] += d0;
			weights[1] -= d0;
		}

		if (knot[dim] + 2 < size[dim]) {
			weights[3] += 0.5f * d1;
			weights[1] -= 0.5f * d1;
		} else {
			weights[2] += d1;
			weights[1] -= d1;
		}
	}

	Float result = 0.0f;
	for (int z=-1; z<=2; ++z) {
		Float wz = knotWeights[2][z+1];
		for (int y=-1; y<=2; ++y) {
			Float wyz = knotWeights[1][y+1] * wz;
			for (int x=-1; x<=2; ++x) {
				Float wxyz = knotWeights[0][x+1] * wyz;

				if (wxyz == 0)
					continue;

				size_t pos = ((knot[2] + z) * size[1] + (knot[1] + y))
					* size[0] + knot[0] + x;

				result += data[pos] * wxyz;
			}
		}
	}
	return result;
}

Float interpCubic3DIrregular(const Point3 &p, const Float **nodes_,
			const Float *data, const Size3 &size, bool extrapolate) {
	Float knotWeights[3][4];
	Size3 knot;

	/* Compute interpolation weights separately for each dimension */
	for (int dim=0; dim<3; ++dim) {
		const Float *nodes = nodes_[dim];
		Float *weights = knotWeights[dim];

		/* Give up when given an out-of-range or NaN argument */
		if (!(p[dim] >= nodes[0] && p[dim] <= nodes[size[dim]-1]) && !extrapolate)
			return 0.0f;

		/* Find the index of the left knot in the queried subinterval, be
		   robust to cases where 't' lies exactly on the right endpoint */
		size_t k = (size_t) std::max((ptrdiff_t) 0, std::min((ptrdiff_t) size[dim] - 2,
			std::lower_bound(nodes, nodes + size[dim], p[dim]) - nodes - 1));
		knot[dim] = k;

		Float width = nodes[k+1] - nodes[k], invWidth = 1 / width;

		/* Compute the relative position within the interval */
		Float t = (p[dim] - nodes[k]) * invWidth,
			  t2 = t*t, t3 = t2*t;

		/* Compute node weights */
		weights[0] = 0.0f;
		weights[1] = 2*t3 - 3*t2 + 1;
		weights[2] = -2*t3 + 3*t2;
		weights[3] = 0.0f;

		/* Derivative weights */
		Float d0 = (t3 - 2*t2 + t) * width,
			  d1 = (t3 - t2) * width;

		/* Turn derivative weights into node weights using
		   an appropriate chosen finite differences stencil */
		if (k > 0) {
			Float factor = 1 / (nodes[k+1]-nodes[k-1]);
			weights[2] += d0 * factor;
			weights[0] -= d0 * factor;
		} else {
			weights[2] += d0 * invWidth;
			weights[1] -= d0 * invWidth;
		}

		if (k + 2 < size[dim]) {
			Float factor = 1 / (nodes[k+2]-nodes[k]);
			weights[3] += d1 * factor;
			weights[1] -= d1 * factor;
		} else {
			weights[2] += d1 * invWidth;
			weights[1] -= d1 * invWidth;
		}
	}

	Float result = 0.0f;
	for (int z=-1; z<=2; ++z) {
		Float wz = knotWeights[2][z+1];
		for (int y=-1; y<=2; ++y) {
			Float wyz = knotWeights[1][y+1] * wz;
			for (int x=-1; x<=2; ++x) {
				Float wxyz = knotWeights[0][x+1] * wyz;

				if (wxyz == 0)
					continue;

				size_t pos = ((knot[2] + z) * size[1] + (knot[1] + y))
					* size[0] + knot[0] + x;

				result += data[pos] * wxyz;
			}
		}
	}
	return result;
}

void stratifiedSample1D(Random *random, Float *dest, int count, bool jitter) {
	Float invCount = 1.0f / count;

	for (int i=0; i<count; i++) {
		Float offset = jitter ? random->nextFloat() : 0.5f;
		*dest++ = (i + offset) * invCount;
	}
}

void stratifiedSample2D(Random *random, Point2 *dest, int countX, int countY, bool jitter) {
	Float invCountX = 1.0f / countX;
	Float invCountY = 1.0f / countY;

	for (int x=0; x<countX; x++) {
		for (int y=0; y<countY; y++) {
			Float offsetX = jitter ? random->nextFloat() : 0.5f;
			Float offsetY = jitter ? random->nextFloat() : 0.5f;
			*dest++ = Point2(
				(x + offsetX) * invCountX,
				(y + offsetY) * invCountY
			);
		}
	}
}

void latinHypercube(Random *random, Float *dest, size_t nSamples, size_t nDim) {
	Float delta = 1 / (Float) nSamples;
	for (size_t i = 0; i < nSamples; ++i)
		for (size_t j = 0; j < nDim; ++j)
			dest[nDim * i + j] = (i + random->nextFloat()) * delta;
	for (size_t i = 0; i < nDim; ++i) {
		for (size_t j = 0; j < nSamples; ++j) {
			size_t other = random->nextSize(nSamples);
			std::swap(dest[nDim * j + i], dest[nDim * other + i]);
		}
	}
}

Vector sphericalDirection(Float theta, Float phi) {
	Float sinTheta, cosTheta, sinPhi, cosPhi;

	math::sincos(theta, &sinTheta, &cosTheta);
	math::sincos(phi, &sinPhi, &cosPhi);

	return Vector(
		sinTheta * cosPhi,
		sinTheta * sinPhi,
		cosTheta
	);
}

void coordinateSystem(const Vector &a, Vector &b, Vector &c) {
	if (std::abs(a.x) > std::abs(a.y)) {
		Float invLen = 1.0f / std::sqrt(a.x * a.x + a.z * a.z);
		c = Vector(a.z * invLen, 0.0f, -a.x * invLen);
	} else {
		Float invLen = 1.0f / std::sqrt(a.y * a.y + a.z * a.z);
		c = Vector(0.0f, a.z * invLen, -a.y * invLen);
	}
	b = cross(c, a);
}

Point2 toSphericalCoordinates(const Vector &v) {
	Point2 result(
		std::acos(v.z),
		std::atan2(v.y, v.x)
	);
	if (result.y < 0)
		result.y += 2*M_PI;
	return result;
}

Float fresnelDielectric(Float cosThetaI, Float cosThetaT, Float eta) {
	if (EXPECT_NOT_TAKEN(eta == 1))
		return 0.0f;

	Float Rs = (cosThetaI - eta * cosThetaT)
			 / (cosThetaI + eta * cosThetaT);
	Float Rp = (eta * cosThetaI - cosThetaT)
			 / (eta * cosThetaI + cosThetaT);

	/* No polarization -- return the unpolarized reflectance */
	return 0.5f * (Rs * Rs + Rp * Rp);
}

Float fresnelDielectricExt(Float cosThetaI_, Float &cosThetaT_, Float eta) {
	if (EXPECT_NOT_TAKEN(eta == 1)) {
		cosThetaT_ = -cosThetaI_;
		return 0.0f;
	}

	/* Using Snell's law, calculate the squared sine of the
	   angle between the normal and the transmitted ray */
	Float scale = (cosThetaI_ > 0) ? 1/eta : eta,
	      cosThetaTSqr = 1 - (1-cosThetaI_*cosThetaI_) * (scale*scale);

	/* Check for total internal reflection */
	if (cosThetaTSqr <= 0.0f) {
		cosThetaT_ = 0.0f;
		return 1.0f;
	}

	/* Find the absolute cosines of the incident/transmitted rays */
	Float cosThetaI = std::abs(cosThetaI_);
	Float cosThetaT = std::sqrt(cosThetaTSqr);

	Float Rs = (cosThetaI - eta * cosThetaT)
			 / (cosThetaI + eta * cosThetaT);
	Float Rp = (eta * cosThetaI - cosThetaT)
			 / (eta * cosThetaI + cosThetaT);

	cosThetaT_ = (cosThetaI_ > 0) ? -cosThetaT : cosThetaT;

	/* No polarization -- return the unpolarized reflectance */
	return 0.5f * (Rs * Rs + Rp * Rp);
}

Spectrum fresnelConductor(Float cosThetaI, const Spectrum &eta, const Spectrum &k) {
	Spectrum tmp = (eta*eta + k*k) * (cosThetaI * cosThetaI);

	Spectrum rParl2 = (tmp - (eta * (2.0f * cosThetaI)) + Spectrum(1.0f))
					/ (tmp + (eta * (2.0f * cosThetaI)) + Spectrum(1.0f));

	Spectrum tmpF = eta*eta + k*k;

	Spectrum rPerp2 = (tmpF - (eta * (2.0f * cosThetaI)) + Spectrum(cosThetaI*cosThetaI)) /
					  (tmpF + (eta * (2.0f * cosThetaI)) + Spectrum(cosThetaI*cosThetaI));

	return (rParl2 + rPerp2) / 2.0f;
}

Vector reflect(const Vector &wi, const Normal &n) {
	return 2 * dot(wi, n) * Vector(n) - wi;
}

Vector refract(const Vector &wi, const Normal &n, Float eta, Float cosThetaT) {
	if (cosThetaT < 0)
		eta = 1 / eta;

	return n * (dot(wi, n) * eta + cosThetaT) - wi * eta;
}

Vector refract(const Vector &wi, const Normal &n, Float eta) {
	if (EXPECT_NOT_TAKEN(eta == 1))
		return -wi;

	Float cosThetaI = dot(wi, n);
	if (cosThetaI > 0)
		eta = 1 / eta;

	/* Using Snell's law, calculate the squared sine of the
	   angle between the normal and the transmitted ray */
	Float cosThetaTSqr = 1 - (1-cosThetaI*cosThetaI) * (eta*eta);

	/* Check for total internal reflection */
	if (cosThetaTSqr <= 0.0f)
		return Vector(0.0f);

	return n * (cosThetaI * eta - math::signum(cosThetaI)
		* std::sqrt(cosThetaTSqr)) - wi * eta;
}

Vector refract(const Vector &wi, const Normal &n, Float eta, Float &cosThetaT, Float &F) {
	Float cosThetaI = dot(wi, n);
	F = fresnelDielectricExt(cosThetaI, cosThetaT, eta);

	if (F == 1.0f) /* Total internal reflection */
		return Vector(0.0f);

	if (cosThetaT < 0)
		eta = 1 / eta;

	return n * (eta * cosThetaI + cosThetaT) - wi * eta;
}

namespace {
	/// Integrand used by fresnelDiffuseReflectance
	inline Float fresnelDiffuseIntegrand(Float eta, Float xi) {
		return fresnelDielectricExt(std::sqrt(xi), eta);
	}
};

Float fresnelDiffuseReflectance(Float eta, bool fast) {
	if (fast) {
		/* Fast mode: the following code approximates the
		 * diffuse Frensel reflectance for the eta<1 and
		 * eta>1 cases. An evalution of the accuracy led
		 * to the following scheme, which cherry-picks
		 * fits from two papers where they are best.
		 */
		if (eta < 1) {
			/* Fit by Egan and Hilgeman (1973). Works
			   reasonably well for "normal" IOR values (<2).

			   Max rel. error in 1.0 - 1.5 : 0.1%
			   Max rel. error in 1.5 - 2   : 0.6%
			   Max rel. error in 2.0 - 5   : 9.5%
			*/
			return -1.4399f * (eta * eta)
				  + 0.7099f * eta
				  + 0.6681f
				  + 0.0636f / eta;
		} else {
			/* Fit by d'Eon and Irving (2011)
			 *
			 * Maintains a good accuracy even for
			 * unrealistic IOR values.
			 *
			 * Max rel. error in 1.0 - 2.0   : 0.1%
			 * Max rel. error in 2.0 - 10.0  : 0.2%
			 */
			Float invEta = 1.0f / eta,
				  invEta2 = invEta*invEta,
				  invEta3 = invEta2*invEta,
				  invEta4 = invEta3*invEta,
				  invEta5 = invEta4*invEta;

			return 0.919317f - 3.4793f * invEta
				 + 6.75335f * invEta2
				 - 7.80989f * invEta3
				 + 4.98554f * invEta4
				 - 1.36881f * invEta5;
		}
	} else {
		GaussLobattoIntegrator quad(1024, 0, 1e-5f);
		return quad.integrate(
			boost::bind(&fresnelDiffuseIntegrand, eta, _1), 0, 1);
	}

	return 0.0f;
}

Float fresnelDiffuseReflectanceSecondMoment(Float eta) {
	/* Fit by d'Eon and Irving (2011) */
	Float invEta = 1.0f / eta,
			  invEta2 = invEta*invEta,
			  invEta3 = invEta2*invEta,
			  invEta4 = invEta3*invEta,
			  invEta5 = invEta4*invEta;
	if (eta < 1) {
		return -1641.1f + 1213.67f * invEta
			 - 568.556f * invEta2
			 + 164.798f * invEta3
			 - 27.0181f * invEta4
			 + 1.91826f * invEta5
			 + 135.926 * eta * eta * eta
			 - 656.175 * eta * eta
			 + 1376.53 * eta;
	} else {
		return 0.828421f - 2.62051f * invEta
			 + 3.36231f * invEta2
			 - 1.95284f * invEta3
			 + 0.236494f * invEta4
			 - 0.145787f * invEta5;
	}
}

std::string timeString(Float time, bool precise) {
	std::ostringstream os;

	if (std::isnan(time) || std::isinf(time))
		return "inf";

	char suffix = 's';
	if (time > 60) {
		time /= 60; suffix = 'm';
		if (time > 60) {
			time /= 60; suffix = 'h';
			if (time > 12) {
				time /= 12; suffix = 'd';
			}
		}
	}

	os << std::setprecision(precise ? 4 : 1)
	   << std::fixed << time << suffix;

	return os.str();
}

Float hypot2(Float a, Float b) {
	Float r;
	if (std::abs(a) > std::abs(b)) {
		r = b / a;
		r = std::abs(a) * std::sqrt(1 + r*r);
	} else if (b != 0) {
		r = a / b;
		r = std::abs(b) * std::sqrt(1 + r*r);
	} else {
		r = 0;
	}
	return r;
}

MTS_NAMESPACE_END
