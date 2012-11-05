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

#pragma once
#if !defined(__MITSUBA_CORE_UTIL_H_)
#define __MITSUBA_CORE_UTIL_H_

#include <boost/static_assert.hpp>

MTS_NAMESPACE_BEGIN

/*! \addtogroup libcore */
/*! @{ */

// -----------------------------------------------------------------------
//! @{ \name String-related utility functions
// -----------------------------------------------------------------------

/**
 * \brief Given a list of delimiters, tokenize
 * a std::string into a vector of strings
 */
extern MTS_EXPORT_CORE std::vector<std::string> tokenize(
	const std::string &string,
	const std::string &delim
);

/// Trim spaces (' ', '\\n', '\\r', '\\t') from the ends of a string
extern MTS_EXPORT_CORE std::string trim(const std::string& str);

/// Indent a string (Used for recursive toString() structure dumping)
extern MTS_EXPORT_CORE std::string indent(const std::string &string, int amount=1);

/// Wrapped snprintf
extern MTS_EXPORT_CORE std::string formatString(const char *pFmt, ...);

/**
 * \brief Convert a time difference (in ms) to a string representation
 * \param time Time value in milliseconds
 * \param precise When set to true, a higher-precision string representation
 * is generated.
 */
extern MTS_EXPORT_CORE std::string timeString(Float time, bool precise = false);

/// Turn a memory size into a human-readable string
extern MTS_EXPORT_CORE std::string memString(size_t size);

/// Return a string representation of a list of objects
template<class Iterator> std::string containerToString(const Iterator &start, const Iterator &end) {
	std::ostringstream oss;
	oss << "{" << std::endl;
	Iterator it = start;
	while (it != end) {
		oss << "  " << indent((*it)->toString());
		++it;
		if (it != end)
			oss << "," << std::endl;
		else
			oss << std::endl;
	}
	oss << "}";
	return oss.str();
}

/// Simple functor for sorting string parameters by length and content
struct SimpleStringOrdering {
	bool operator()(const std::string &a, const std::string &b) const {
		if (a.length() == b.length())
			return a < b;
		return a.length() < b.length();
	}
};

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Miscellaneous
// -----------------------------------------------------------------------

/// Allocate an aligned region of memory
extern MTS_EXPORT_CORE void * __restrict allocAligned(size_t size);

/// Free an aligned region of memory
extern MTS_EXPORT_CORE void freeAligned(void *ptr);

#if defined(WIN32)
/// Return a string version of GetLastError()
extern std::string MTS_EXPORT_CORE lastErrorText();
#endif

/// Determine the number of available CPU cores
extern MTS_EXPORT_CORE int getCoreCount();

/// Return the host name of this machine
extern MTS_EXPORT_CORE std::string getHostName();

/// Return the fully qualified domain name of this machine
extern MTS_EXPORT_CORE std::string getFQDN();

/**
 * \brief Enable floating point exceptions (to catch NaNs, overflows,
 * arithmetic with infinity).
 *
 * On Intel processors, this applies to both x87 and SSE2 math
 *
 * \return \c true if floating point exceptions were active
 * before calling the function
 */
extern MTS_EXPORT_CORE bool enableFPExceptions();

/**
 * \brief Disable floating point exceptions
 *
 * \return \c true if floating point exceptions were active
 * before calling the function
 */
extern MTS_EXPORT_CORE bool disableFPExceptions();

/// Restore floating point exceptions to the specified state
extern MTS_EXPORT_CORE void restoreFPExceptions(bool state);

/// Cast between types that have an identical binary representation.
template<typename T, typename U> inline T union_cast(const U &val) {
	BOOST_STATIC_ASSERT(sizeof(T) == sizeof(U));

	union {
		U u;
		T t;
	} caster = {val};

	return caster.t;
}

/// Swaps the byte order of the underlying representation
template<typename T> inline T endianness_swap(T value) {
	union {
		T value;
		uint8_t byteValue[sizeof(T)];
	} u;

	u.value = value;
	std::reverse(&u.byteValue[0], &u.byteValue[sizeof(T)]);
	return u.value;
}

#ifdef __GNUC__
#if defined(__i386__)
static FINLINE uint64_t rdtsc(void) {
  uint64_t x;
	 __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
	 return x;
}
#elif defined(__x86_64__)
static FINLINE uint64_t rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t) lo)| (((uint64_t) hi) << 32);
}
#endif
#elif defined(__MSVC__)
static FINLINE __int64 rdtsc(void) {
	return __rdtsc();
}
#else
# error "Cannot generate the rdtsc intrinsic."
#endif

/**
 * \brief Apply an arbitrary permutation to an array in linear time
 *
 * This algorithm is based on Donald Knuth's book
 * "The Art of Computer Programming, Volume 3: Sorting and Searching"
 * (1st edition, section 5.2, page 595)
 *
 * Given a permutation and an array of values, it applies the permutation
 * in linear time without requiring additional memory. This is based on
 * the fact that each permutation can be decomposed into a disjoint set
 * of permutations, which can then be applied individually.
 *
 * \param data
 *     Pointer to the data that should be permuted
 * \param perm
 *     Input permutation vector having the same size as \c data. After
 *     the function terminates, this vector will be set to the
 *     identity permutation.
 */
template <typename DataType, typename IndexType> void permute_inplace(
		DataType *data, std::vector<IndexType> &perm) {
	for (size_t i=0; i<perm.size(); i++) {
		if (perm[i] != i) {
			/* The start of a new cycle has been found. Save
			   the value at this position, since it will be
			   overwritten */
			IndexType j = (IndexType) i;
			DataType curval = data[i];

			do {
				/* Shuffle backwards */
				IndexType k = perm[j];
				data[j] = data[k];

				/* Also fix the permutations on the way */
				perm[j] = j;
				j = k;

				/* Until the end of the cycle has been found */
			} while (perm[j] != i);

			/* Fix the final position with the saved value */
			data[j] = curval;
			perm[j] = j;
		}
	}
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Numerical utility functions
// -----------------------------------------------------------------------

/// sqrt(a^2 + b^2) without underflow (like 'hypot' on compilers that support C99)
extern MTS_EXPORT_CORE Float hypot2(Float a, Float b);

/// Base-2 logarithm
extern MTS_EXPORT_CORE Float log2(Float value);

/// Always-positive modulo function (assumes b > 0)
inline int modulo(int a, int b) {
	int r = a % b;
	return (r < 0) ? r+b : r;
}

/// Always-positive modulo function, float version (assumes b > 0)
inline Float modulo(Float a, Float b) {
	Float r = std::fmod(a, b);
	return (r < 0) ? r+b : r;
}

/// Compute the signum (a.k.a. "sign") function
inline Float signum(Float value) {
	if (value < 0)
		return -1;
	else if (value > 0)
		return 1;
	else return 0;
}

/// Integer floor function
inline int floorToInt(Float value) {
	return (int) std::floor(value);
}

/// Integer ceil function
inline int ceilToInt(Float value) {
	return (int) std::ceil(value);
}

/// Base-2 logarithm (32-bit integer version)
extern MTS_EXPORT_CORE int log2i(uint32_t value);

/// Base-2 logarithm (64-bit integer version)
extern MTS_EXPORT_CORE int log2i(uint64_t value);

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline int log2i(size_t value) {
	if (sizeof(size_t) == 8)
		return log2i((uint64_t) value);
	else
		return log2i((uint32_t) value);
}
#endif

/// Check if an integer is a power of two (unsigned 32 bit version)
inline bool isPowerOfTwo(uint32_t i) { return (i & (i-1)) == 0; }

/// Check if an integer is a power of two (signed 32 bit version)
inline bool isPowerOfTwo(int32_t i) { return i > 0 && (i & (i-1)) == 0; }

/// Check if an integer is a power of two (64 bit version)
inline bool isPowerOfTwo(uint64_t i) { return (i & (i-1)) == 0; }

/// Check if an integer is a power of two (signed 64 bit version)
inline bool isPowerOfTwo(int64_t i) { return i > 0 && (i & (i-1)) == 0; }

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline bool isPowerOfTwo(size_t value) {
	if (sizeof(size_t) == 8) /// will be optimized away
		return isPowerOfTwo((uint64_t) value);
	else
		return isPowerOfTwo((uint32_t) value);
}
#endif

/// Round an integer to the next power of two
extern MTS_EXPORT_CORE uint32_t roundToPowerOfTwo(uint32_t i);

/// Round an integer to the next power of two (64 bit version)
extern MTS_EXPORT_CORE uint64_t roundToPowerOfTwo(uint64_t i);

#if defined(MTS_AMBIGUOUS_SIZE_T)
/// Round an integer to the next power of two
inline size_t roundToPowerOfTwo(size_t value) {
	if (sizeof(size_t) == 8) /// will be optimized away
		return (size_t) roundToPowerOfTwo((uint64_t) value);
	else
		return (size_t) roundToPowerOfTwo((uint32_t) value);
}
#endif

/**
 * \brief Solve a quadratic equation of the form a*x^2 + b*x + c = 0.
 * \return \c true if a solution could be found
 */
extern MTS_EXPORT_CORE bool solveQuadratic(Float a, Float b,
	Float c, Float &x0, Float &x1);

/**
 * \brief Solve a double-precision quadratic equation of the
 * form a*x^2 + b*x + c = 0.
 * \return \c true if a solution could be found
 */
extern MTS_EXPORT_CORE bool solveQuadraticDouble(double a, double b,
	double c, double &x0, double &x1);

/**
 * \brief Evaluate a cubic spline interpolant of a regularly sampled 1D function
 *
 * This implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param x
 *      Evaluation point
 * \param data
 *      Floating point array containing \c size regularly spaced evaluations
 *      in the range [\c min,\c max] of the function to be approximated.
 * \param min
 *      Position of the first knot
 * \param max
 *      Position of the last knot
 * \param size
 *      Denotes the size of the \c data array
 * \param extrapolate
 *      Extrapolate data values when \c x is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt>
 *      and \c x lies outside of [\c min, \c max]
 */
extern MTS_EXPORT_CORE Float interpCubic1D(Float x, const Float *data,
		Float min, Float max, size_t size, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of an \a irregularly sampled 1D function
 *
 * This implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param x
 *      Evaluation point
 * \param nodes
 *      Floating point array containing \c size irregularly spaced values
 *      denoting positions the where the function to be interpolated was evaluated.
 *      They must be provided in \a increasing order.
 * \param data
 *      Floating point array containing interpolant values matched to
 *      the entries of \c nodes.
 * \param size
 *      Denotes the size of the \c data array
 * \param extrapolate
 *      Extrapolate data values when \c x is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt>
 *      and \c x lies outside of \a [\c min, \c max]
 */
extern MTS_EXPORT Float interpCubic1DIrregular(Float x, const Float *nodes,
		const Float *data, size_t size, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of a regularly sampled 2D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * patch.
 *
 * \param p
 *      Evaluation point
 * \param data
 *      A 2D floating point array of <tt>size.x*size.y</tt> cells containing regularly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate.
 * \param min
 *      Position of the first knot on each dimension
 * \param max
 *      Position of the last knot on each dimension
 * \param size
 *      Denotes the size of the \c data array (along each dimension)
 * \param extrapolate
 *      Extrapolate data values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float interpCubic2D(const Point2 &p, const Float *data,
		const Point2 &min, const Point2 &max, const Size2 &size, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of an \a irregularly sampled 2D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * region.
 *
 * When the underlying function is sampled on a regular grid, \ref interpCubic2D()
 * should be preferred, since data lookups will be considerably faster.
 *
 * \param p
 *      Evaluation point
 * \param nodes
 *      Pointer to a list for each dimension denoting the positions where the function
 *      to be interpolated was evaluated. The <tt>i</tt>-th array must have
 *      size <tt>size[i]</tt> and contain position values in \a increasing order.
 * \param data
 *      A 2D floating point array of <tt>size.x*size.y</tt> cells containing irregularly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate.
 * \param size
 *      Denotes the size of the \c data array (along each dimension)
 * \param extrapolate
 *      Extrapolate data values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float interpCubic2DIrregular(const Point2 &p, const Float **nodes,
		const Float *data, const Size2 &size, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of a regularly sampled 3D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * region.
 *
 * \param p
 *      Evaluation point of the interpolant
 * \param data
 *      A 3D floating point array of <tt>size.x*size.y*size.z</tt> cells containing regularly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate,
 *      then 'y', and finally 'z' increments.
 * \param min
 *      Position of the first knot on each dimension
 * \param max
 *      Position of the last knot on each dimension
 * \param size
 *      Denotes the size of the \c data array (along each dimension)
 * \param extrapolate
 *      Extrapolate data values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float interpCubic3D(const Point3 &p, const Float *data,
		const Point3 &min, const Point3 &max, const Size3 &size, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of an \a irregularly sampled 3D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * region.
 *
 * When the underlying function is sampled on a regular grid, \ref interpCubic3D()
 * should be preferred, since data lookups will be considerably faster.
 *
 * \param p
 *      Evaluation point
 * \param nodes
 *      Pointer to a list for each dimension denoting the positions where the function
 *      to be interpolated was evaluated. The <tt>i</tt>-th array must have
 *      size <tt>size[i]</tt> and contain position values in \a increasing order.
 * \param data
 *      A 2D floating point array of <tt>size.x*size.y</tt> cells containing irregularly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate,
 *      then 'y', and finally 'z' increments.
 * \param size
 *      Denotes the size of the \c data array (along each dimension)
 * \param extrapolate
 *      Extrapolate data values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float interpCubic3DIrregular(const Point3 &p, const Float **nodes,
		const Float *data, const Size3 &size, bool extrapolate = false);

//// Convert radians to degrees
inline Float radToDeg(Float value) { return value * (180.0f / M_PI); }

/// Convert degrees to radians
inline Float degToRad(Float value) { return value * (M_PI / 180.0f); }

/// Generic clamping function
template <typename Scalar> inline Scalar clamp(Scalar value, Scalar min, Scalar max) {
	return std::min(max, std::max(min, value));
}

/// Linearly interpolate between two values
inline Float lerp(Float t, Float v1, Float v2) {
    return ((Float) 1 - t) * v1 + t * v2;
}

/// S-shaped smoothly varying interpolation between two values
inline Float smoothStep(Float min, Float max, Float value) {
    Float v = clamp((value - min) / (max - min), (Float) 0, (Float) 1);
    return v * v * (-2 * v  + 3);
}

/**
 * \brief Numerically well-behaved routine for computing the angle
 * between two unit direction vectors
 *
 * This should be used wherever one is tempted to compute the
 * arc cosine of a dot product!
 *
 * Proposed by Don Hatch at
 * http://www.plunk.org/~hatch/rightway.php
 */
template <typename VectorType> inline Float unitAngle(const VectorType &u, const VectorType &v) {
	if (dot(u, v) < 0)
		return M_PI - 2 * std::asin(0.5f * (v+u).length());
	else
		return 2 * std::asin(0.5f * (v-u).length());
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Warping and sampling-related utility functions
// -----------------------------------------------------------------------

/**
 * \brief Solve a 2x2 linear equation system using basic linear algebra
 */
extern MTS_EXPORT_CORE bool solveLinearSystem2x2(const Float a[2][2], const Float b[2], Float x[2]);

/**
 * \brief Complete the set {a} to an orthonormal base
 * \remark In Python, this function is used as
 *     follows: <tt>s, t = coordinateSystem(n)</tt>
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE void coordinateSystem(const Vector &a, Vector &b, Vector &c);

/**
 * \brief Generate (optionally jittered) stratified 1D samples
 * \param random Source of random numbers
 * \param dest A pointer to a floating point array with at least
 *             count entries
 * \param count The interval [0, 1] is split into count strata
 * \param jitter Randomly jitter the samples?
 */
extern MTS_EXPORT_CORE void stratifiedSample1D(Random *random, Float *dest,
	int count, bool jitter);

/**
 * \brief Generate (optionally jittered) stratified 2D samples
 * \param random Source of random numbers
 * \param dest A pointer to a floating point array
 * \param countX The X axis interval [0, 1] is split into countX strata
 * \param countY The Y axis interval [0, 1] is split into countY strata
 * \param jitter Randomly jitter the samples?
 */
extern MTS_EXPORT_CORE void stratifiedSample2D(Random *random, Point2 *dest,
	int countX, int countY, bool jitter);

/// Generate latin hypercube samples
extern MTS_EXPORT_CORE void latinHypercube(
		Random *random, Float *dest, size_t nSamples, size_t nDim);

/// Convert spherical coordinates to a direction
extern MTS_EXPORT_CORE Vector sphericalDirection(Float theta, Float phi);

/// Convert a direction to spherical coordinates
extern MTS_EXPORT_CORE Point2 toSphericalCoordinates(const Vector &v);

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Fresnel reflectance computation and related things
// -----------------------------------------------------------------------

/**
 * \brief Calculates the unpolarized Fresnel reflection coefficient
 * at a planar interface between two dielectrics
 *
 * This is a basic implementation that just returns the value of
 * \f[
 * R(\cos\theta_i,\cos\theta_t,\eta)=\frac{1}{2} \left[
 * \left(\frac{\eta\cos\theta_i-\cos\theta_t}{\eta\cos\theta_i+\cos\theta_t}\right)^2+
 * \left(\frac{\cos\theta_i-\eta\cos\theta_t}{\cos\theta_i+\eta\cos\theta_t}\right)^2
 * \right]
 * \f]
 * The transmitted direction must be provided. There is no logic pertaining to
 * total internal reflection or negative direction cosines.
 *
 * \param cosThetaI
 * 		Absolute cosine of the angle between the normal and the incident ray
 * \param cosThetaT
 * 		Absolute cosine of the angle between the normal and the transmitted ray
 * \param eta
 * 		Relative refractive index to the transmitted direction
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Float fresnelDielectric(Float cosThetaI,
		Float cosThetaT, Float eta);

/**
 * \brief Calculates the unpolarized Fresnel reflection coefficient
 * at a planar interface between two dielectrics (extended version)
 *
 * In comparison to \ref fresnelDielectric(), this function internally
 * computes the transmitted direction and returns it using the \c cosThetaT
 * argument. When encountering total internal reflection, it sets
 * <tt>cosThetaT=0</tt> and returns the value 1.
 *
 * When <tt>cosThetaI < 0</tt>, the function computes the Fresnel reflectance
 * from the \a internal boundary, which is equivalent to calling the function
 * with arguments <tt>fresnelDielectric(abs(cosThetaI), cosThetaT, 1/eta)</tt>.
 *
 * \remark When accessed from Python, this function has the signature
 * "<tt>F, cosThetaT = fresnelDielectricExt(cosThetaI, eta)</tt>".
 *
 * \param cosThetaI
 * 		Cosine of the angle between the normal and the incident ray
 * 		(may be negative)
 * \param cosThetaT
 * 		Argument used to return the cosine of the angle between the normal
 * 		and the transmitted ray, will have the opposite sign of \c cosThetaI
 * \param eta
 * 		Relative refractive index
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Float fresnelDielectricExt(Float cosThetaI,
	Float &cosThetaT, Float eta);

/**
 * \brief Calculates the unpolarized Fresnel reflection coefficient
 * at a planar interface between two dielectrics (extended version)
 *
 * This is just a convenience wrapper function around the other \c fresnelDielectricExt
 * function, which does not return the transmitted direction cosine in case it is
 * not needed by the application.
 *
 * \param cosThetaI
 * 		Cosine of the angle between the normal and the incident ray
 */
inline Float fresnelDielectricExt(Float cosThetaI, Float eta) { Float cosThetaT;
	return fresnelDielectricExt(cosThetaI, cosThetaT, eta); }

/**
 * \brief Calculates the unpolarized fresnel reflection coefficient
 * at a planar interface between vacuum and a conductor.
 *
 * \param cosThetaI
 * 		Cosine of the angle between the normal and the incident ray
 * \param eta
 * 		Real refractive index (wavelength-dependent)
 * \param k
 * 		Imaginary refractive index (wavelength-dependent)
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Spectrum fresnelConductor(Float cosThetaI,
		const Spectrum &eta, const Spectrum &k);

/**
 * \brief Calculates the diffuse unpolarized fresnel reflectance of
 * a dielectric material (sometimes referred to as "Fdr").
 *
 * This value quantifies what fraction of diffuse incident illumination
 * will, on average, be reflected at a dielectric material boundary
 *
 * \param eta
 *      Relative refraction coefficient
 * \param fast
 *      Compute an approximate value? If set to \c true, the
 *      implementation will use a polynomial approximation with
 *      a max relative error of ~0.5% on the interval 0.5 < \c eta < 2.
 *      When \c fast=false, the code will use Gauss-Lobatto quadrature
 *      to compute the diffuse reflectance more accurately, and for
 *      a wider range of refraction coefficients, but at a cost
 *      in terms of performance.
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Float fresnelDiffuseReflectance(
	Float eta, bool fast = false);

/**
 * \brief Computes the second moment of the Fresnel reflectance
 * of a dielectric material
 *
 * This function computes the integral between the Fresnel reflectance
 * and the squared cosine of the angle of incidence. The function
 * \ref fresnelDiffuseReflectance() correspondingly computes what
 * can be thought of as the first moment.
 */
extern MTS_EXPORT_CORE Float fresnelDiffuseReflectanceSecondMoment(Float eta);

/**
 * \brief Specularly reflect direction \c wi with respect to the given surface normal
 * \param wi
 *     Incident direction
 * \param n
 *     Surface normal
 * \return
 *     Specularly reflected direction
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Vector reflect(const Vector &wi, const Normal &n);

/**
 * \brief Specularly refract the direction \c wi into a planar dielectric with
 * the given surface normal and index of refraction.
 *
 * This variant internally computes the transmitted direction cosine by
 * calling \ref fresnelDielectricExt. As a side result, the cosine and
 * Fresnel reflectance are computed and returned via the reference arguments
 * \c cosThetaT and \c F.
 *
 * \remark When accessed from Python, this function has the signature
 * "<tt>dir, cosThetaT, F = refract(wi, n, eta)</tt>".
 *
 * \param wi
 *     Incident direction
 * \param n
 *     Surface normal
 * \param eta
 *     Relative index of refraction at the interface
 * \param cosThetaT
 *     Parameter used to return the signed cosine of the angle between the transmitted
 *     direction and the surface normal
 * \param F
 *     Parameter used to return the Fresnel reflectance
 * \return
 *     Specularly transmitted direction (or zero in
 *     the case of total internal reflection)
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Vector refract(const Vector &wi, const Normal &n,
	Float eta, Float &cosThetaT, Float &F);

/**
 * \brief Specularly refract the direction \c wi into a planar dielectric with
 * the given surface normal and index of refraction.
 *
 * This variant assumes that the transmitted direction cosine has
 * has <em>already</em> been computed, allowing it to save some time.
 *
 * \param wi
 *     Incident direction
 * \param n
 *     Surface normal
 * \param eta
 *     Relative index of refraction at the interface
 * \param cosThetaT
 *     Signed cosine of the angle between the transmitted direction and
 *     the surface normal obtained from a prior call to \ref fresnelDielectricExt()
 * \return
 *     Specularly transmitted direction
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Vector refract(const Vector &wi, const Normal &n,
	Float eta, Float cosThetaT);

/**
 * \brief Specularly refract the direction \c wi into a planar dielectric with
 * the given surface normal and index of refraction.
 *
 * This function is a simple convenience function that only returns the refracted
 * direction while not computing the Frensel reflectance.
 *
 * \param wi
 *     Incident direction
 * \param n
 *     Surface normal
 * \param eta
 *     Relative index of refraction at the interface
 * \return
 *     Specularly transmitted direction (or zero in
 *     the case of total internal reflection)
 * \ingroup libpython
 */
extern MTS_EXPORT_CORE Vector refract(const Vector &wi, const Normal &n, Float eta);

//! @}
// -----------------------------------------------------------------------

/*! @} */

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_UTIL_H_ */
