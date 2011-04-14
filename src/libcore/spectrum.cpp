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

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

void Spectrum::staticInitialization() {
	if (SPECTRUM_SAMPLES == 3) {
		/* CIE RGB primaries */
		m_wavelengths[0] = 700.0f;
		m_wavelengths[1] = 546.1f;
		m_wavelengths[2] = 438.8f;
	} else {
		SAssert(SPECTRUM_SAMPLES <= CIE_count);
		std::ostringstream oss;
		Float stepSize = SPECTRUM_RANGE / (Float) SPECTRUM_SAMPLES;
		for (int i=0; i<SPECTRUM_SAMPLES; i++) {
			int wavelength = SPECTRUM_MIN_WAVELENGTH
							+ (int) (stepSize * i);
			m_wavelengths[i] = (Float) wavelength;
			oss << wavelength;
			if (i+1<SPECTRUM_SAMPLES)
				oss << ", ";
		}
		Spectrum s(1.0f);
		m_normalization = s.getLuminance();
		m_invNormalization = 1.0f / m_normalization;
		SLog(EDebug, "Configuring for spectral rendering using "
			"the wavelengths %s", oss.str().c_str());
	}
}

void Spectrum::staticShutdown() {
	/* Do nothing */
}

void Spectrum::fromSmoothSpectrum(const SmoothSpectrum *smooth) {
#if SPECTRUM_SAMPLES == 3
	int index = 0;
	Float x=0, y=0, z=0, ySum=0;
	for (int lambda=CIE_start; lambda<=CIE_end; ++lambda) {
		const Float value = smooth->eval((Float) lambda);
		x += CIE_X[index] * value;
		y += CIE_Y[index] * value;
		z += CIE_Z[index] * value;
		ySum += CIE_Y[index];
		++index;
	}
	x /= ySum; y /= ySum; z /= ySum;
	fromXYZ(x, y, z);
#else
	for (int i=0; i<SPECTRUM_SAMPLES; i++)
		s[i] = smooth->eval(m_wavelengths[i]);
#endif
}

Float Spectrum::eval(Float lambda) const {
#if SPECTRUM_SAMPLES == 3
	SLog(EError, "Spectrum::eval() is not supported when "
		"using RGB rendering");
#else
	const Float pos = SPECTRUM_SAMPLES * (lambda - 
		SPECTRUM_MIN_WAVELENGTH) / (Float) SPECTRUM_RANGE;
	int index = (int) pos;
	const Float offset = pos - index;
	if (index < 0 || index >= SPECTRUM_SAMPLES)
		return 0.0f;

	if (index == SPECTRUM_SAMPLES-1)
		return s[index];
	else 
		return (1.0f - offset) * s[index] + offset * s[index+1];
#endif
	return 0.0f;
}

#if SPECTRUM_SAMPLES == 3
void Spectrum::fromXYZ(Float x, Float y, Float z) {
	s[0] =  3.240479f * x + -1.537150f * y + -0.498535f * z;
	s[1] = -0.969256f * x +  1.875991f * y +  0.041556f * z;
	s[2] =  0.055648f * x + -0.204043f * y +  1.057311f * z;
}

void Spectrum::toXYZ(Float &x, Float &y, Float &z) const {
	x = s[0] * 0.412453f + s[1] * 0.357580f + s[2] * 0.180423f;
	y = s[0] * 0.212671f + s[1] * 0.715160f + s[2] * 0.072169f;
	z = s[0] * 0.019334f + s[1] * 0.119193f + s[2] * 0.950227f;
}

#else
void Spectrum::toXYZ(Float &x, Float &y, Float &z) const {
	int index = 0;
	x = y = z = 0.0f;
	for (int lambda=CIE_start; lambda<=CIE_end; ++lambda) {
		const Float value = eval(lambda);
		x += CIE_X[index] * value;
		y += CIE_Y[index] * value;
		z += CIE_Z[index] * value;
		++index;
	}
	x *= CIE_normalization;
	y *= CIE_normalization;
	z *= CIE_normalization;
}

Float Spectrum::getLuminance() const {
	Float luminance = 0.0f;
	int index = 0;
	for (int lambda=CIE_start; lambda<=CIE_end; ++lambda) 
		luminance += (CIE_Y[index++] * eval(lambda)) * CIE_normalization;
	return luminance;
}

void Spectrum::fromXYZ(Float x, Float y, Float z) {
	/* Convert XYZ to linear RGB according to ITU-R Rec. BT.709 */
	Float r = 3.240479f*x -1.537150f*y - 0.498535f*z;
	Float g = -0.969256f*x + 1.875991f*y + 0.041556f*z;
	Float b = 0.055648f*x - 0.204043f*y + 1.057311f*z;

	/* Method suggested by Shirley and Marschner */
	for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
		Float wavelength = m_wavelengths[i];
		if (wavelength < 400)
			s[i] = 0.0f;
		else if (wavelength < 500)
			s[i] = b;
		else if (wavelength < 600)
			s[i] = g;
		else if (wavelength <= 700)
			s[i] = r;
		else
			s[i] = 0.0f;
	}
}

#endif

#if SPECTRUM_SAMPLES != 3
void Spectrum::toLinearRGB(Float &r, Float &g, Float &b) const {
	Float x, y, z;
	toXYZ(x, y, z);
	/* Convert XYZ to linear RGB according to ITU-R Rec. BT.709 */
	r = std::max((Float) 0.0f, m_invNormalization * (3.240479f*x -1.537150f*y - 0.498535f*z));
	g = std::max((Float) 0.0f, m_invNormalization * (-0.969256f*x + 1.875991f*y + 0.041556f*z));
	b = std::max((Float) 0.0f, m_invNormalization * (0.055648f*x - 0.204043f*y + 1.057311f*z));
}

void Spectrum::fromLinearRGB(Float r, Float g, Float b) {
	/* Convert linear RGB to XYZ according to ITU-R Rec. BT.709 */
	Float x = 0.412453f*r + 0.35758f*g + 0.180423f*b;
	Float y = 0.212671f*r + 0.71516f*g + 0.0721688f*b;
	Float z = 0.0193338f*r + 0.119194f*g + 0.950227f*b;
	fromXYZ(x, y, z);
}
#endif

inline Float toSRGBComponent(Float value) {
	if (value <= (Float) 0.0031308)
		return (Float) 12.92 * value;
	return (Float) (1.0 + 0.055)
		* std::pow(value, (Float) (1.0/2.4))
		 - (Float) 0.055;
}

inline Float fromSRGBComponent(Float value) {
	if (value <= (Float) 0.04045)
		return value / (Float) 12.92;
	return std::pow((value + (Float) 0.055)
		/ (Float) (1.0 + 0.055), (Float) 2.4);
}

void Spectrum::toSRGB(Float &r, Float &g, Float &b) const {
	toLinearRGB(r, g, b);
	r = toSRGBComponent(r);
	g = toSRGBComponent(g);
	b = toSRGBComponent(b);
}

void Spectrum::fromSRGB(Float r, Float g, Float b) {
	r = fromSRGBComponent(r);
	g = fromSRGBComponent(g);
	b = fromSRGBComponent(b);
	fromLinearRGB(r, g, b);
}

/* RGBE conversion based on Bruce Walter's and Greg Ward's code
  -> http://www.graphics.cornell.edu/online/formats/rgbe/rgbe.c */
void Spectrum::toRGBE(uint8_t rgbe[4]) const {
	Float r, g, b;
	toLinearRGB(r, g, b);
	/* Find the largest contribution */
	Float max = std::max(std::max(r, g), b);
	if (max < 1e-32) {
		rgbe[0] = rgbe[1] = rgbe[2] = rgbe[3] = 0;
	} else {
		int e;
		/* Extract exponent and convert the fractional part into
		 the [0..255] range. Afterwards, divide by max so that
		 any color component multiplied by the result will be in [0,255] */
		max = std::frexp(max, &e) * (Float) 256 / max;
		rgbe[0] = (uint8_t) (r * max);
		rgbe[1] = (uint8_t) (g * max);
		rgbe[2] = (uint8_t) (b * max);
		rgbe[3] = e+128; /* Exponent value in bias format */
	}
}

/* RGBE conversion based on Bruce Walter's and Greg Ward's code
  -> http://www.graphics.cornell.edu/online/formats/rgbe/rgbe.c */
void Spectrum::fromRGBE(const uint8_t rgbe[4]) {
	if (rgbe[3]) {
		/* Calculate exponent/256 */
		Float exp = std::ldexp((Float) 1, (int) rgbe[3] - (128+8));
		fromLinearRGB(rgbe[0]*exp, rgbe[1]*exp, rgbe[2]*exp);
	} else {
		s[0] = s[1] = s[2] = 0.0f;
	}
}

std::string Spectrum::toString() const {
	std::ostringstream oss;
	oss << "[";
	for (int i=0; i<SPECTRUM_SAMPLES; i++) {
#if SPECTRUM_SAMPLES == 3
		oss << s[i];
#else
		oss << m_wavelengths[i] << ":" << s[i];
#endif
		if (i < SPECTRUM_SAMPLES - 1)
			oss << ", ";
	}
	oss << "]";
	return oss.str();
}

Float BlackBodySpectrum::eval(Float l) const {
	/* Convert inputs to meters and kelvins */
	const double lambda = l * 1e-9;
	const double c = 299792458;
	const double h = 6.62607095e-34;
	const double k = 1.3807e-23;
	/* Watts per unit surface area (cm^-2) per unit wavelength (nm^-1) */
	const double I = (2*M_PI*h*c*c) * std::pow(lambda, -5.0)
		/ (1e9*1e4*(std::exp((h/k)*c/(lambda*m_temperature)) - 1.0));
	return (Float) I;
}

void InterpolatedSpectrum::appendSample(Float lambda, Float value) {
	SAssert(m_wavelength.size() == 0 || m_wavelength[m_wavelength.size()-1] < lambda);
	m_wavelength.push_back(lambda);
	m_value.push_back(value);
}

Float InterpolatedSpectrum::eval(Float lambda) const {
	typedef std::vector<Float>::const_iterator it;
	SAssert(m_wavelength.size() > 0);
	if (lambda < m_wavelength[0] || lambda > m_wavelength[m_wavelength.size()-1])
		return 0;
	std::pair<it, it> result = 
		std::equal_range(m_wavelength.begin(), m_wavelength.end(), lambda);

	size_t idx1 = result.first - m_wavelength.begin();
	size_t idx2 = result.second - m_wavelength.begin();
	if (idx1 == idx2) {
		Float x_a = m_wavelength[idx1-1];
		Float x_b = m_wavelength[idx1];
		Float y_a = m_value[idx1-1];
		Float y_b = m_value[idx1];
		Float t = (lambda - x_a) / (x_b-x_a);
		return y_a + t * (y_b-y_a);
	} else if (idx2 == idx1+1) {
		/* Hit a value exactly */
		return m_value[idx1];
	} else {
		SLog(EError, "Internal error while interpolating spectrum values");
		return 0;
	}
}

Float Spectrum::m_wavelengths[SPECTRUM_SAMPLES];
Float Spectrum::m_normalization;
Float Spectrum::m_invNormalization;
const Float Spectrum::CIE_normalization = 683.002f;

/* CIE 1931 XYZ color matching functions from 
   http://www.cvrl.org/database/data/cmfs/ciexyz31_1.txt */
const Float Spectrum::CIE_X[Spectrum::CIE_count] = {
	0.0001299000f, 0.0001458470f, 0.0001638021f, 0.0001840037f,
	0.0002066902f, 0.0002321000f, 0.0002607280f, 0.0002930750f,
	0.0003293880f, 0.0003699140f, 0.0004149000f, 0.0004641587f,
	0.0005189860f, 0.0005818540f, 0.0006552347f, 0.0007416000f,
	0.0008450296f, 0.0009645268f, 0.001094949f, 0.001231154f,
	0.001368000f, 0.001502050f, 0.001642328f, 0.001802382f,
	0.001995757f, 0.002236000f, 0.002535385f, 0.002892603f,
	0.003300829f, 0.003753236f, 0.004243000f, 0.004762389f,
	0.005330048f, 0.005978712f, 0.006741117f, 0.007650000f,
	0.008751373f, 0.01002888f, 0.01142170f, 0.01286901f,
	0.01431000f, 0.01570443f, 0.01714744f, 0.01878122f,
	0.02074801f, 0.02319000f, 0.02620736f, 0.02978248f,
	0.03388092f, 0.03846824f, 0.04351000f, 0.04899560f,
	0.05502260f, 0.06171880f, 0.06921200f, 0.07763000f,
	0.08695811f, 0.09717672f, 0.1084063f, 0.1207672f,
	0.1343800f, 0.1493582f, 0.1653957f, 0.1819831f,
	0.1986110f, 0.2147700f, 0.2301868f, 0.2448797f,
	0.2587773f, 0.2718079f, 0.2839000f, 0.2949438f,
	0.3048965f, 0.3137873f, 0.3216454f, 0.3285000f,
	0.3343513f, 0.3392101f, 0.3431213f, 0.3461296f,
	0.3482800f, 0.3495999f, 0.3501474f, 0.3500130f,
	0.3492870f, 0.3480600f, 0.3463733f, 0.3442624f,
	0.3418088f, 0.3390941f, 0.3362000f, 0.3331977f,
	0.3300411f, 0.3266357f, 0.3228868f, 0.3187000f,
	0.3140251f, 0.3088840f, 0.3032904f, 0.2972579f,
	0.2908000f, 0.2839701f, 0.2767214f, 0.2689178f,
	0.2604227f, 0.2511000f, 0.2408475f, 0.2298512f,
	0.2184072f, 0.2068115f, 0.1953600f, 0.1842136f,
	0.1733273f, 0.1626881f, 0.1522833f, 0.1421000f,
	0.1321786f, 0.1225696f, 0.1132752f, 0.1042979f,
	0.09564000f, 0.08729955f, 0.07930804f, 0.07171776f,
	0.06458099f, 0.05795001f, 0.05186211f, 0.04628152f,
	0.04115088f, 0.03641283f, 0.03201000f, 0.02791720f,
	0.02414440f, 0.02068700f, 0.01754040f, 0.01470000f,
	0.01216179f, 0.009919960f, 0.007967240f, 0.006296346f,
	0.004900000f, 0.003777173f, 0.002945320f, 0.002424880f,
	0.002236293f, 0.002400000f, 0.002925520f, 0.003836560f,
	0.005174840f, 0.006982080f, 0.009300000f, 0.01214949f,
	0.01553588f, 0.01947752f, 0.02399277f, 0.02910000f,
	0.03481485f, 0.04112016f, 0.04798504f, 0.05537861f,
	0.06327000f, 0.07163501f, 0.08046224f, 0.08973996f,
	0.09945645f, 0.1096000f, 0.1201674f, 0.1311145f,
	0.1423679f, 0.1538542f, 0.1655000f, 0.1772571f,
	0.1891400f, 0.2011694f, 0.2133658f, 0.2257499f,
	0.2383209f, 0.2510668f, 0.2639922f, 0.2771017f,
	0.2904000f, 0.3038912f, 0.3175726f, 0.3314384f,
	0.3454828f, 0.3597000f, 0.3740839f, 0.3886396f,
	0.4033784f, 0.4183115f, 0.4334499f, 0.4487953f,
	0.4643360f, 0.4800640f, 0.4959713f, 0.5120501f,
	0.5282959f, 0.5446916f, 0.5612094f, 0.5778215f,
	0.5945000f, 0.6112209f, 0.6279758f, 0.6447602f,
	0.6615697f, 0.6784000f, 0.6952392f, 0.7120586f,
	0.7288284f, 0.7455188f, 0.7621000f, 0.7785432f,
	0.7948256f, 0.8109264f, 0.8268248f, 0.8425000f,
	0.8579325f, 0.8730816f, 0.8878944f, 0.9023181f,
	0.9163000f, 0.9297995f, 0.9427984f, 0.9552776f,
	0.9672179f, 0.9786000f, 0.9893856f, 0.9995488f,
	1.0090892f, 1.0180064f, 1.0263000f, 1.0339827f,
	1.0409860f, 1.0471880f, 1.0524667f, 1.0567000f,
	1.0597944f, 1.0617992f, 1.0628068f, 1.0629096f,
	1.0622000f, 1.0607352f, 1.0584436f, 1.0552244f,
	1.0509768f, 1.0456000f, 1.0390369f, 1.0313608f,
	1.0226662f, 1.0130477f, 1.0026000f, 0.9913675f,
	0.9793314f, 0.9664916f, 0.9528479f, 0.9384000f,
	0.9231940f, 0.9072440f, 0.8905020f, 0.8729200f,
	0.8544499f, 0.8350840f, 0.8149460f, 0.7941860f,
	0.7729540f, 0.7514000f, 0.7295836f, 0.7075888f,
	0.6856022f, 0.6638104f, 0.6424000f, 0.6215149f,
	0.6011138f, 0.5811052f, 0.5613977f, 0.5419000f,
	0.5225995f, 0.5035464f, 0.4847436f, 0.4661939f,
	0.4479000f, 0.4298613f, 0.4120980f, 0.3946440f,
	0.3775333f, 0.3608000f, 0.3444563f, 0.3285168f,
	0.3130192f, 0.2980011f, 0.2835000f, 0.2695448f,
	0.2561184f, 0.2431896f, 0.2307272f, 0.2187000f,
	0.2070971f, 0.1959232f, 0.1851708f, 0.1748323f,
	0.1649000f, 0.1553667f, 0.1462300f, 0.1374900f,
	0.1291467f, 0.1212000f, 0.1136397f, 0.1064650f,
	0.09969044f, 0.09333061f, 0.08740000f, 0.08190096f,
	0.07680428f, 0.07207712f, 0.06768664f, 0.06360000f,
	0.05980685f, 0.05628216f, 0.05297104f, 0.04981861f,
	0.04677000f, 0.04378405f, 0.04087536f, 0.03807264f,
	0.03540461f, 0.03290000f, 0.03056419f, 0.02838056f,
	0.02634484f, 0.02445275f, 0.02270000f, 0.02108429f,
	0.01959988f, 0.01823732f, 0.01698717f, 0.01584000f,
	0.01479064f, 0.01383132f, 0.01294868f, 0.01212920f,
	0.01135916f, 0.01062935f, 0.009938846f, 0.009288422f,
	0.008678854f, 0.008110916f, 0.007582388f, 0.007088746f,
	0.006627313f, 0.006195408f, 0.005790346f, 0.005409826f,
	0.005052583f, 0.004717512f, 0.004403507f, 0.004109457f,
	0.003833913f, 0.003575748f, 0.003334342f, 0.003109075f,
	0.002899327f, 0.002704348f, 0.002523020f, 0.002354168f,
	0.002196616f, 0.002049190f, 0.001910960f, 0.001781438f,
	0.001660110f, 0.001546459f, 0.001439971f, 0.001340042f,
	0.001246275f, 0.001158471f, 0.001076430f, 0.0009999493f,
	0.0009287358f, 0.0008624332f, 0.0008007503f, 0.0007433960f,
	0.0006900786f, 0.0006405156f, 0.0005945021f, 0.0005518646f,
	0.0005124290f, 0.0004760213f, 0.0004424536f, 0.0004115117f,
	0.0003829814f, 0.0003566491f, 0.0003323011f, 0.0003097586f,
	0.0002888871f, 0.0002695394f, 0.0002515682f, 0.0002348261f,
	0.0002191710f, 0.0002045258f, 0.0001908405f, 0.0001780654f,
	0.0001661505f, 0.0001550236f, 0.0001446219f, 0.0001349098f,
	0.0001258520f, 0.0001174130f, 0.0001095515f, 0.0001022245f,
	0.00009539445f, 0.00008902390f, 0.00008307527f, 0.00007751269f,
	0.00007231304f, 0.00006745778f, 0.00006292844f, 0.00005870652f,
	0.00005477028f, 0.00005109918f, 0.00004767654f, 0.00004448567f,
	0.00004150994f, 0.00003873324f, 0.00003614203f, 0.00003372352f,
	0.00003146487f, 0.00002935326f, 0.00002737573f, 0.00002552433f,
	0.00002379376f, 0.00002217870f, 0.00002067383f, 0.00001927226f,
	0.00001796640f, 0.00001674991f, 0.00001561648f, 0.00001455977f,
	0.00001357387f, 0.00001265436f, 0.00001179723f, 0.00001099844f,
	0.00001025398f, 0.000009559646f, 0.000008912044f, 0.000008308358f,
	0.000007745769f, 0.000007221456f, 0.000006732475f, 0.000006276423f,
	0.000005851304f, 0.000005455118f, 0.000005085868f, 0.000004741466f,
	0.000004420236f, 0.000004120783f, 0.000003841716f, 0.000003581652f,
	0.000003339127f, 0.000003112949f, 0.000002902121f, 0.000002705645f,
	0.000002522525f, 0.000002351726f, 0.000002192415f, 0.000002043902f,
	0.000001905497f, 0.000001776509f, 0.000001656215f, 0.000001544022f,
	0.000001439440f, 0.000001341977f, 0.000001251141f
};

const Float Spectrum::CIE_Y[Spectrum::CIE_count] = {
	0.000003917000f, 0.000004393581f, 0.000004929604f, 0.000005532136f,
	0.000006208245f, 0.000006965000f, 0.000007813219f, 0.000008767336f,
	0.000009839844f, 0.00001104323f, 0.00001239000f, 0.00001388641f,
	0.00001555728f, 0.00001744296f, 0.00001958375f, 0.00002202000f,
	0.00002483965f, 0.00002804126f, 0.00003153104f, 0.00003521521f,
	0.00003900000f, 0.00004282640f, 0.00004691460f, 0.00005158960f,
	0.00005717640f, 0.00006400000f, 0.00007234421f, 0.00008221224f,
	0.00009350816f, 0.0001061361f, 0.0001200000f, 0.0001349840f,
	0.0001514920f, 0.0001702080f, 0.0001918160f, 0.0002170000f,
	0.0002469067f, 0.0002812400f, 0.0003185200f, 0.0003572667f,
	0.0003960000f, 0.0004337147f, 0.0004730240f, 0.0005178760f,
	0.0005722187f, 0.0006400000f, 0.0007245600f, 0.0008255000f,
	0.0009411600f, 0.001069880f, 0.001210000f, 0.001362091f,
	0.001530752f, 0.001720368f, 0.001935323f, 0.002180000f,
	0.002454800f, 0.002764000f, 0.003117800f, 0.003526400f,
	0.004000000f, 0.004546240f, 0.005159320f, 0.005829280f,
	0.006546160f, 0.007300000f, 0.008086507f, 0.008908720f,
	0.009767680f, 0.01066443f, 0.01160000f, 0.01257317f,
	0.01358272f, 0.01462968f, 0.01571509f, 0.01684000f,
	0.01800736f, 0.01921448f, 0.02045392f, 0.02171824f,
	0.02300000f, 0.02429461f, 0.02561024f, 0.02695857f,
	0.02835125f, 0.02980000f, 0.03131083f, 0.03288368f,
	0.03452112f, 0.03622571f, 0.03800000f, 0.03984667f,
	0.04176800f, 0.04376600f, 0.04584267f, 0.04800000f,
	0.05024368f, 0.05257304f, 0.05498056f, 0.05745872f,
	0.06000000f, 0.06260197f, 0.06527752f, 0.06804208f,
	0.07091109f, 0.07390000f, 0.07701600f, 0.08026640f,
	0.08366680f, 0.08723280f, 0.09098000f, 0.09491755f,
	0.09904584f, 0.1033674f, 0.1078846f, 0.1126000f,
	0.1175320f, 0.1226744f, 0.1279928f, 0.1334528f,
	0.1390200f, 0.1446764f, 0.1504693f, 0.1564619f,
	0.1627177f, 0.1693000f, 0.1762431f, 0.1835581f,
	0.1912735f, 0.1994180f, 0.2080200f, 0.2171199f,
	0.2267345f, 0.2368571f, 0.2474812f, 0.2586000f,
	0.2701849f, 0.2822939f, 0.2950505f, 0.3085780f,
	0.3230000f, 0.3384021f, 0.3546858f, 0.3716986f,
	0.3892875f, 0.4073000f, 0.4256299f, 0.4443096f,
	0.4633944f, 0.4829395f, 0.5030000f, 0.5235693f,
	0.5445120f, 0.5656900f, 0.5869653f, 0.6082000f,
	0.6293456f, 0.6503068f, 0.6708752f, 0.6908424f,
	0.7100000f, 0.7281852f, 0.7454636f, 0.7619694f,
	0.7778368f, 0.7932000f, 0.8081104f, 0.8224962f,
	0.8363068f, 0.8494916f, 0.8620000f, 0.8738108f,
	0.8849624f, 0.8954936f, 0.9054432f, 0.9148501f,
	0.9237348f, 0.9320924f, 0.9399226f, 0.9472252f,
	0.9540000f, 0.9602561f, 0.9660074f, 0.9712606f,
	0.9760225f, 0.9803000f, 0.9840924f, 0.9874812f,
	0.9903128f, 0.9928116f, 0.9949501f, 0.9967108f,
	0.9980983f, 0.9991120f, 0.9997482f, 1.0000000f,
	0.9998567f, 0.9993046f, 0.9983255f, 0.9968987f,
	0.9950000f, 0.9926005f, 0.9897426f, 0.9864444f,
	0.9827241f, 0.9786000f, 0.9740837f, 0.9691712f,
	0.9638568f, 0.9581349f, 0.9520000f, 0.9454504f,
	0.9384992f, 0.9311628f, 0.9234576f, 0.9154000f,
	0.9070064f, 0.8982772f, 0.8892048f, 0.8797816f,
	0.8700000f, 0.8598613f, 0.8493920f, 0.8386220f,
	0.8275813f, 0.8163000f, 0.8047947f, 0.7930820f,
	0.7811920f, 0.7691547f, 0.7570000f, 0.7447541f,
	0.7324224f, 0.7200036f, 0.7074965f, 0.6949000f,
	0.6822192f, 0.6694716f, 0.6566744f, 0.6438448f,
	0.6310000f, 0.6181555f, 0.6053144f, 0.5924756f,
	0.5796379f, 0.5668000f, 0.5539611f, 0.5411372f,
	0.5283528f, 0.5156323f, 0.5030000f, 0.4904688f,
	0.4780304f, 0.4656776f, 0.4534032f, 0.4412000f,
	0.4290800f, 0.4170360f, 0.4050320f, 0.3930320f,
	0.3810000f, 0.3689184f, 0.3568272f, 0.3447768f,
	0.3328176f, 0.3210000f, 0.3093381f, 0.2978504f,
	0.2865936f, 0.2756245f, 0.2650000f, 0.2547632f,
	0.2448896f, 0.2353344f, 0.2260528f, 0.2170000f,
	0.2081616f, 0.1995488f, 0.1911552f, 0.1829744f,
	0.1750000f, 0.1672235f, 0.1596464f, 0.1522776f,
	0.1451259f, 0.1382000f, 0.1315003f, 0.1250248f,
	0.1187792f, 0.1127691f, 0.1070000f, 0.1014762f,
	0.09618864f, 0.09112296f, 0.08626485f, 0.08160000f,
	0.07712064f, 0.07282552f, 0.06871008f, 0.06476976f,
	0.06100000f, 0.05739621f, 0.05395504f, 0.05067376f,
	0.04754965f, 0.04458000f, 0.04175872f, 0.03908496f,
	0.03656384f, 0.03420048f, 0.03200000f, 0.02996261f,
	0.02807664f, 0.02632936f, 0.02470805f, 0.02320000f,
	0.02180077f, 0.02050112f, 0.01928108f, 0.01812069f,
	0.01700000f, 0.01590379f, 0.01483718f, 0.01381068f,
	0.01283478f, 0.01192000f, 0.01106831f, 0.01027339f,
	0.009533311f, 0.008846157f, 0.008210000f, 0.007623781f,
	0.007085424f, 0.006591476f, 0.006138485f, 0.005723000f,
	0.005343059f, 0.004995796f, 0.004676404f, 0.004380075f,
	0.004102000f, 0.003838453f, 0.003589099f, 0.003354219f,
	0.003134093f, 0.002929000f, 0.002738139f, 0.002559876f,
	0.002393244f, 0.002237275f, 0.002091000f, 0.001953587f,
	0.001824580f, 0.001703580f, 0.001590187f, 0.001484000f,
	0.001384496f, 0.001291268f, 0.001204092f, 0.001122744f,
	0.001047000f, 0.0009765896f, 0.0009111088f, 0.0008501332f,
	0.0007932384f, 0.0007400000f, 0.0006900827f, 0.0006433100f,
	0.0005994960f, 0.0005584547f, 0.0005200000f, 0.0004839136f,
	0.0004500528f, 0.0004183452f, 0.0003887184f, 0.0003611000f,
	0.0003353835f, 0.0003114404f, 0.0002891656f, 0.0002684539f,
	0.0002492000f, 0.0002313019f, 0.0002146856f, 0.0001992884f,
	0.0001850475f, 0.0001719000f, 0.0001597781f, 0.0001486044f,
	0.0001383016f, 0.0001287925f, 0.0001200000f, 0.0001118595f,
	0.0001043224f, 0.00009733560f, 0.00009084587f, 0.00008480000f,
	0.00007914667f, 0.00007385800f, 0.00006891600f, 0.00006430267f,
	0.00006000000f, 0.00005598187f, 0.00005222560f, 0.00004871840f,
	0.00004544747f, 0.00004240000f, 0.00003956104f, 0.00003691512f,
	0.00003444868f, 0.00003214816f, 0.00003000000f, 0.00002799125f,
	0.00002611356f, 0.00002436024f, 0.00002272461f, 0.00002120000f,
	0.00001977855f, 0.00001845285f, 0.00001721687f, 0.00001606459f,
	0.00001499000f, 0.00001398728f, 0.00001305155f, 0.00001217818f,
	0.00001136254f, 0.00001060000f, 0.000009885877f, 0.000009217304f,
	0.000008592362f, 0.000008009133f, 0.000007465700f, 0.000006959567f,
	0.000006487995f, 0.000006048699f, 0.000005639396f, 0.000005257800f,
	0.000004901771f, 0.000004569720f, 0.000004260194f, 0.000003971739f,
	0.000003702900f, 0.000003452163f, 0.000003218302f, 0.000003000300f,
	0.000002797139f, 0.000002607800f, 0.000002431220f, 0.000002266531f,
	0.000002113013f, 0.000001969943f, 0.000001836600f, 0.000001712230f,
	0.000001596228f, 0.000001488090f, 0.000001387314f, 0.000001293400f,
	0.000001205820f, 0.000001124143f, 0.000001048009f, 0.0000009770578f,
	0.0000009109300f, 0.0000008492513f, 0.0000007917212f, 0.0000007380904f,
	0.0000006881098f, 0.0000006415300f, 0.0000005980895f, 0.0000005575746f,
	0.0000005198080f, 0.0000004846123f, 0.0000004518100f
};
const Float Spectrum::CIE_Z[Spectrum::CIE_count] = {
	0.0006061000f, 0.0006808792f, 0.0007651456f, 0.0008600124f,
	0.0009665928f, 0.001086000f, 0.001220586f, 0.001372729f,
	0.001543579f, 0.001734286f, 0.001946000f, 0.002177777f,
	0.002435809f, 0.002731953f, 0.003078064f, 0.003486000f,
	0.003975227f, 0.004540880f, 0.005158320f, 0.005802907f,
	0.006450001f, 0.007083216f, 0.007745488f, 0.008501152f,
	0.009414544f, 0.01054999f, 0.01196580f, 0.01365587f,
	0.01558805f, 0.01773015f, 0.02005001f, 0.02251136f,
	0.02520288f, 0.02827972f, 0.03189704f, 0.03621000f,
	0.04143771f, 0.04750372f, 0.05411988f, 0.06099803f,
	0.06785001f, 0.07448632f, 0.08136156f, 0.08915364f,
	0.09854048f, 0.1102000f, 0.1246133f, 0.1417017f,
	0.1613035f, 0.1832568f, 0.2074000f, 0.2336921f,
	0.2626114f, 0.2947746f, 0.3307985f, 0.3713000f,
	0.4162091f, 0.4654642f, 0.5196948f, 0.5795303f,
	0.6456000f, 0.7184838f, 0.7967133f, 0.8778459f,
	0.9594390f, 1.0390501f, 1.1153673f, 1.1884971f,
	1.2581233f, 1.3239296f, 1.3856000f, 1.4426352f,
	1.4948035f, 1.5421903f, 1.5848807f, 1.6229600f,
	1.6564048f, 1.6852959f, 1.7098745f, 1.7303821f,
	1.7470600f, 1.7600446f, 1.7696233f, 1.7762637f,
	1.7804334f, 1.7826000f, 1.7829682f, 1.7816998f,
	1.7791982f, 1.7758671f, 1.7721100f, 1.7682589f,
	1.7640390f, 1.7589438f, 1.7524663f, 1.7441000f,
	1.7335595f, 1.7208581f, 1.7059369f, 1.6887372f,
	1.6692000f, 1.6475287f, 1.6234127f, 1.5960223f,
	1.5645280f, 1.5281000f, 1.4861114f, 1.4395215f,
	1.3898799f, 1.3387362f, 1.2876400f, 1.2374223f,
	1.1878243f, 1.1387611f, 1.0901480f, 1.0419000f,
	0.9941976f, 0.9473473f, 0.9014531f, 0.8566193f,
	0.8129501f, 0.7705173f, 0.7294448f, 0.6899136f,
	0.6521049f, 0.6162000f, 0.5823286f, 0.5504162f,
	0.5203376f, 0.4919673f, 0.4651800f, 0.4399246f,
	0.4161836f, 0.3938822f, 0.3729459f, 0.3533000f,
	0.3348578f, 0.3175521f, 0.3013375f, 0.2861686f,
	0.2720000f, 0.2588171f, 0.2464838f, 0.2347718f,
	0.2234533f, 0.2123000f, 0.2011692f, 0.1901196f,
	0.1792254f, 0.1685608f, 0.1582000f, 0.1481383f,
	0.1383758f, 0.1289942f, 0.1200751f, 0.1117000f,
	0.1039048f, 0.09666748f, 0.08998272f, 0.08384531f,
	0.07824999f, 0.07320899f, 0.06867816f, 0.06456784f,
	0.06078835f, 0.05725001f, 0.05390435f, 0.05074664f,
	0.04775276f, 0.04489859f, 0.04216000f, 0.03950728f,
	0.03693564f, 0.03445836f, 0.03208872f, 0.02984000f,
	0.02771181f, 0.02569444f, 0.02378716f, 0.02198925f,
	0.02030000f, 0.01871805f, 0.01724036f, 0.01586364f,
	0.01458461f, 0.01340000f, 0.01230723f, 0.01130188f,
	0.01037792f, 0.009529306f, 0.008749999f, 0.008035200f,
	0.007381600f, 0.006785400f, 0.006242800f, 0.005749999f,
	0.005303600f, 0.004899800f, 0.004534200f, 0.004202400f,
	0.003900000f, 0.003623200f, 0.003370600f, 0.003141400f,
	0.002934800f, 0.002749999f, 0.002585200f, 0.002438600f,
	0.002309400f, 0.002196800f, 0.002100000f, 0.002017733f,
	0.001948200f, 0.001889800f, 0.001840933f, 0.001800000f,
	0.001766267f, 0.001737800f, 0.001711200f, 0.001683067f,
	0.001650001f, 0.001610133f, 0.001564400f, 0.001513600f,
	0.001458533f, 0.001400000f, 0.001336667f, 0.001270000f,
	0.001205000f, 0.001146667f, 0.001100000f, 0.001068800f,
	0.001049400f, 0.001035600f, 0.001021200f, 0.001000000f,
	0.0009686400f, 0.0009299200f, 0.0008868800f, 0.0008425600f,
	0.0008000000f, 0.0007609600f, 0.0007236800f, 0.0006859200f,
	0.0006454400f, 0.0006000000f, 0.0005478667f, 0.0004916000f,
	0.0004354000f, 0.0003834667f, 0.0003400000f, 0.0003072533f,
	0.0002831600f, 0.0002654400f, 0.0002518133f, 0.0002400000f,
	0.0002295467f, 0.0002206400f, 0.0002119600f, 0.0002021867f,
	0.0001900000f, 0.0001742133f, 0.0001556400f, 0.0001359600f,
	0.0001168533f, 0.0001000000f, 0.00008613333f, 0.00007460000f,
	0.00006500000f, 0.00005693333f, 0.00004999999f, 0.00004416000f,
	0.00003948000f, 0.00003572000f, 0.00003264000f, 0.00003000000f,
	0.00002765333f, 0.00002556000f, 0.00002364000f, 0.00002181333f,
	0.00002000000f, 0.00001813333f, 0.00001620000f, 0.00001420000f,
	0.00001213333f, 0.00001000000f, 0.000007733333f, 0.000005400000f,
	0.000003200000f, 0.000001333333f, 0.000000000000f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f
};

MTS_NAMESPACE_END
