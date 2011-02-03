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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Kajiya-Kay phase function
 *
 * The function is normalized so that it has no energy loss when 
 * \a ks=1 and illumination arrives perpendicularly to the surface.
 */
class KajiyaKayPhaseFunction : public PhaseFunction {
public:
	KajiyaKayPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
		m_ks = props.getFloat("ks", .4f);
		m_kd = props.getFloat("kd", .2f);
		m_exponent = props.getFloat("exponent", 4.0f);
		if (m_kd + m_ks > 1.0f)
			Log(EWarn, "Energy conservation is violated!");
	}

	KajiyaKayPhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
		m_ks = stream->readFloat();
		m_kd = stream->readFloat();
		m_exponent = stream->readFloat();
		configure();
	}

	virtual ~KajiyaKayPhaseFunction() { }

	void configure() {
		/* Compute the normalization for perpendicular illumination
		   using Simpson quadrature */
		int nParts = 1000;
		Float stepSize = M_PI / nParts, m=4, theta = stepSize;

		m_normalization = 0; /* 0 at the endpoints */
		for (int i=1; i<nParts; ++i) {
			Float value = std::pow(std::cos(theta - M_PI/2), m_exponent)
				* std::sin(theta);
			m_normalization += value * m;
			theta += stepSize;
			m = 6-m;
		}

		m_normalization = 1/(m_normalization * stepSize/3 * 2 * M_PI);
		Log(EDebug, "Kajiya-kay normalization factor = %f", m_normalization);
	}


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeFloat(m_ks);
		stream->writeFloat(m_kd);
		stream->writeFloat(m_exponent);
	}

	Spectrum sample(PhaseFunctionQueryRecord &pRec,
			Sampler *sampler) const {
		pRec.wo = squareToSphere(sampler->next2D());
		return f(pRec) * (4 * M_PI);
	}

	Spectrum sample(PhaseFunctionQueryRecord &pRec, 
			Float &pdf, Sampler *sampler) const {
		pRec.wo = squareToSphere(sampler->next2D());
		pdf = 1/(4 * M_PI);
		return f(pRec);
	}

	Float pdf(const PhaseFunctionQueryRecord &pRec) const {
		return 1/(4 * M_PI);
	}

	Spectrum f(const PhaseFunctionQueryRecord &pRec) const {
		if (pRec.mRec.orientation.length() == 0)
			return Spectrum(m_kd / (4*M_PI));

		Frame frame(normalize(pRec.mRec.orientation));
		Vector reflectedLocal = frame.toLocal(pRec.wo);

		reflectedLocal.z = -dot(pRec.wi, frame.n);
		Float a = std::sqrt((1-reflectedLocal.z*reflectedLocal.z) / 
			(reflectedLocal.x*reflectedLocal.x + reflectedLocal.y*reflectedLocal.y));
		reflectedLocal.y *= a;
		reflectedLocal.x *= a;
		Vector R = frame.toWorld(reflectedLocal);

		return Spectrum((std::pow(std::max((Float) 0, dot(R, pRec.wo)), m_exponent))
			* m_normalization * m_ks + m_kd / (4*M_PI));
	}

	std::string toString() const {
		return "KajiyaKayPhaseFunction[]";
	}

#if 0
	/// For testing purposes
	Float integrateOverOutgoing(const Vector &wi) {
		MediumSamplingRecord mRec;
		mRec.orientation = Vector(0,0,1);
		int res = 100;
		/* Nested composite Simpson's rule */
		Float hExt = M_PI / res,
		      hInt = (2*M_PI)/(res*2);

		Float result = 0;

		for (int i=0; i<=res; ++i) {
			Float theta = hExt*i;
			Float weightExt = (i & 1) ? 4.0f : 2.0f;
			if (i == 0 || i == res)
				weightExt = 1;

			for (int j=0; j<=res*2; ++j) {
				Float phi = hInt*j;
				Float weightInt = (j & 1) ? 4.0f : 2.0f;
				if (j == 0 || j == 2*res)
					weightInt = 1;

				Float value = f(mRec, wi, sphericalDirection(theta, phi))[0]*std::sin(theta)
					* weightExt*weightInt;

				result += value;
			}
		}
		
		return hExt*hInt/9;
	}
#endif

	MTS_DECLARE_CLASS()
private:
	Float m_ks, m_kd, m_exponent, m_normalization;
};


MTS_IMPLEMENT_CLASS_S(KajiyaKayPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(KajiyaKayPhaseFunction, "Kajiya-Kay phase function");
MTS_NAMESPACE_END
