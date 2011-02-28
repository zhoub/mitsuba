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

#if !defined(__SPECTRUM_H)
#define __SPECTRUM_H

#include <mitsuba/mitsuba.h>

#define SPECTRUM_MIN_WAVELENGTH   400
#define SPECTRUM_MAX_WAVELENGTH   700
#define SPECTRUM_RANGE            (SPECTRUM_MAX_WAVELENGTH-SPECTRUM_MIN_WAVELENGTH+1)
#define SPECTRUM_SAMPLES          3

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract smooth spectral power distribution data type,
 * which supports evaluation at arbitrary wavelengths
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE SmoothSpectrum {
public:
	/**
	 * Evaluate the value of the spectral power distribution
	 * at the given wavelength.
	 */
	virtual Float eval(Float lambda) const = 0;
	
	virtual ~SmoothSpectrum() { }
};

/**
 * \brief Spectral power distribution based on Planck's black body law
 *
 * Computes the spectral power distribution of a black body of the 
 * specified temperature.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE BlackBodySpectrum : public SmoothSpectrum {
public:
	/**
	 * \brief Construct a new black body spectrum given the emitter's
	 * temperature in Kelvin.
	 */
	inline BlackBodySpectrum(Float temperature) {
		m_temperature = temperature;
	}

	virtual ~BlackBodySpectrum() { }

	/** \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 */
	virtual Float eval(Float lambda) const;
private:
	Float m_temperature;
};

/**
 * \brief Linearly interpolated spectral power distribution
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE InterpolatedSpectrum : public SmoothSpectrum {
public:
	/**
	 * \brief Create a new interpolated spectrum with space 
	 * for the specified number of samples
	 */
	inline InterpolatedSpectrum(size_t size = 0) {
		m_wavelength.reserve(size);
		m_value.reserve(size);
	}

	/**
	 * \brief Append an entry to the spectral power distribution.
	 *
	 * Entries must be added in order of increasing wavelength
	 */
	void appendSample(Float lambda, Float value);

	/**
	 * \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 */
	virtual Float eval(Float lambda) const;
	
	virtual ~InterpolatedSpectrum() { }
private:
	std::vector<Float> m_wavelength, m_value;
};

/** \brief Discrete spectral power distribution 
 * based on a (usually small) number of wavelength samples. 
 *
 * When SPECTRUM_SAMPLES is set to 3 (the default), this class 
 * falls back to linear RGB as its internal representation.
 *
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE Spectrum : public Eigen::Array<Float, SPECTRUM_SAMPLES, 1> {
public:
	typedef Eigen::Array<Float, SPECTRUM_SAMPLES, 1> Base;

	/// Create a new spectral power distribution, but don't initialize the contents
#if !defined(MTS_DEBUG_UNINITIALIZED)
	inline Spectrum() : Base() { }
#else
	inline Spectrum() : Base() {
		setConstant(std::numeric_limits<double>::quiet_NaN()))
	}
#endif

	/// Create a new spectral power distribution with all samples set to the given value
	explicit inline Spectrum(Float v) : Base(Spectrum::Constant(v)) { }

	/// Generic copy constructor
    template<typename OtherDerived> Spectrum(
		const Eigen::ArrayBase<OtherDerived>& other) : Base(other) {}

	/// Unserialize a spectral power distribution from a binary data stream
	explicit inline Spectrum(Stream *stream) : Base(stream) { }

    template<typename OtherDerived> Spectrum &operator=
			(const Eigen::ArrayBase<OtherDerived>& other) {
		this->Base::operator=(other);
		return *this;
    }

	/// Returns \c true if the spectrum contains NaN-valued samples
	inline bool hasNaN() const {
		return isnan().any();
	}

	/// Returns \c true if the spectrum only contains valid (non-NaN, nonnegative) samples
	inline bool isValid() const {
		return positive().all();
	}

	/// Return the average over all wavelengths
	inline Float average() const {
		return sum() / SPECTRUM_SAMPLES;
	}

	/// Clamp negative values
	inline void clampNegative() {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			coeffRef(i) = std::max((Float) 0.0f, coeffRef(i));
	}

	/// Check if this spectrum is zero at all wavelengths
	inline bool isZero() const {
		return (array() == 0).all();
	}

	/**
	 * \brief Evaluate the SPD at an arbitrary wavelength
	 * (uses interpolation)
	 */
	Float eval(Float lambda) const;

	/// Return the luminance in candelas.
#if SPECTRUM_SAMPLES == 3
	inline Float getLuminance() const {
		return coeff(0) * 0.212671f + coeff(1) * 0.715160f + coeff(2) * 0.072169f;
	}
#else
	Float getLuminance() const;
#endif

	/// Convert from a spectral power distribution to XYZ colors
	void toXYZ(Float &x, Float &y, Float &z) const;

	/// Convert from XYZ to a spectral power distribution
	void fromXYZ(Float x, Float y, Float z);

#if SPECTRUM_SAMPLES == 3
	/// Convert to linear RGB
	inline void toLinearRGB(Float &r, Float &g, Float &b) const {
		r = coeff(0);
		g = coeff(1);
		b = coeff(2);
	}

	/// Convert from linear RGB
	inline void fromLinearRGB(Float r, Float g, Float b) {
		coeffRef(0) = r;
		coeffRef(1) = g;
		coeffRef(2) = b;
	}
#else
	/// Convert to linear RGB
	void toLinearRGB(Float &r, Float &g, Float &b) const;

	/// Convert from linear RGB
	void fromLinearRGB(Float r, Float g, Float b);	
#endif

	/// Convert to sRGB
	void toSRGB(Float &r, Float &g, Float &b) const;

	/// Convert from sRGB
	void fromSRGB(Float r, Float g, Float b);

	/// Linear RGBE conversion based on Bruce Walter's and Greg Ward's code
	void fromRGBE(const uint8_t rgbe[4]);

	/// Linear RGBE conversion based on Bruce Walter's and Greg Ward's code 
	void toRGBE(uint8_t rgbe[4]) const;

	/// Initialize with spectral values from a smooth spectrum representation
	void fromSmoothSpectrum(const SmoothSpectrum *smooth);

	/// Return the wavelength corresponding to an index
	inline static Float getWavelength(int index) {
		SAssert(index < SPECTRUM_SAMPLES);
		return m_wavelengths[index];
	}

	/// Return a string representation
	std::string toString() const;

	/** 
	 * Static initialization (should be called once during the
	 * application's initialization phase
	 */
	static void staticInitialization();
	static void staticShutdown();
protected:
	/// Configured wavelengths in nanometers
	static Float m_wavelengths[SPECTRUM_SAMPLES];

	/// Normalization factor for XYZ<->RGB conversion
	static Float m_normalization;

	/// Inverse of \ref m_normalization
	static Float m_invNormalization;

	/**
	 * @{ \name CIE 1931 XYZ color matching functions. 
	 * From http://www.cvrl.org/database/data/cmfs/ciexyz31_1.txt
	 */
	static const int   CIE_start = 360;
	static const int   CIE_end   = 830;
	static const int   CIE_count = CIE_end - CIE_start + 1;
	static const Float CIE_normalization;
	static const Float CIE_X[CIE_count];
	static const Float CIE_Y[CIE_count];
	static const Float CIE_Z[CIE_count];
	/// @}
};

MTS_NAMESPACE_END

#endif /* __SPECTRUM_H */
