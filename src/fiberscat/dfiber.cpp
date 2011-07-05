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

#include <mitsuba/render/fiberscat.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/brent.h>
#include <mitsuba/render/sampler.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

/**
 * Diffuse fiber scattering function
 */
class DiffuseFiber : public FiberScatteringFunction {
public:
	DiffuseFiber(const Properties &props) 
		: FiberScatteringFunction(props) {
		m_reflectance = props.getSpectrum("reflectance", Spectrum(.5f));
	}

	DiffuseFiber(Stream *stream, InstanceManager *manager) 
		: FiberScatteringFunction(stream, manager) {
		m_reflectance = Spectrum(stream);
	}

	virtual ~DiffuseFiber() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		FiberScatteringFunction::serialize(stream, manager);
		m_reflectance.serialize(stream);
	}

	static Float cdfErrT(Float xi, Float t) {
		// the objective function for sin(alpha) sampling, described in sampleScattering below.
		return (2/M_PI * (std::acos(1-t) - std::sqrt(2*t-t*t) * (1 - t)) - xi) / std::sqrt(t);
	}

	Spectrum sample(FiberScatteringRecord &pRec, Float &pdf,
			Sampler *sampler) const {
		/* We want to sample the sphere according to cos(theta) where theta is the angle from the normal plane.
		 * We'll do this by choosing a hemisphere and sampling the hemisphere,
		 * so the pdf to sample is (2/pi^2) cos(theta) = (2/pi^2) sin(alpha) where alpha
		 * is the angle to the tangent vector.
		 * The inverse cdf approach says
		 *    phi = 2pi xi1 and xi2 = 2pi int_0^alpha (2/pi^2) sin x (sin x dx)
		 * The result is
		 *    xi2(al) = 2/pi (al - sin(al) cos(al)).
		 * This needs to be solved numerically.  A good problem transformation, arrived at empirically,
		 * is to use the independent variable t = 1 - cos(al) and to solve
		 *    xi2(t) / sqrt(t) = xi2_given / sqrt(t)
		 * using the starting point t = xi2_given^(2/3).
		 */

		Point2 sample = sampler->next2D();
		Float xi1 = sample.x, xi2 = sample.y;

		// Choose hemisphere
		bool posZ = xi2 < 0.5;
		xi2 = posZ ? 2 * xi2 : 2 - 2 * xi2;

		// initial guess
		Float t0 = pow(xi2, 2./3.);

		// solve for t
		// There is an opportunity to use a derivative-aware solver to save time, if needed.
		BrentSolver brentSolver(100, 1e-6f);
		BrentSolver::Result result = brentSolver.solve(boost::bind(cdfErrT, xi2, _1), 0, 1, t0);
		SAssert(result.success);
		Float t = result.x;

		// compute theta and phi
		Float cos_alpha = (posZ ? 1 : -1) * (1 - t);
		Float phi = 2*M_PI * xi1;

		// construct vector
		Float sin_alpha = std::sqrt(1 - cos_alpha*cos_alpha);
		pRec.wo.x = cos_alpha;
		pRec.wo.y = sin_alpha * std::cos(phi);
		pRec.wo.z = sin_alpha * std::sin(phi);

		// set value and pdf; they are constant
		pdf = this->pdf(pRec);
		return evaluate(pRec);
	}

	Spectrum evaluate(const FiberScatteringRecord &pRec) const {
		return m_reflectance / (M_PI*M_PI);
	}

	Float pdf(const FiberScatteringRecord &pRec) const {
		// PDF of wo | wi is cos(th) / pi^2.
		Float cos_theta = std::sqrt(1 - pRec.wo.x*pRec.wo.x);
		return cos_theta / (M_PI*M_PI);
	}

	std::string toString() const {
		return "DiffuseFiber[]";
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_reflectance;
};


/* Matlab hemisphere sampling code for reference:

function [alpha, phi] = sin_sample(xi1, xi2)

% Given (x1, x2) uniformly distributed in [0,1]^2, produce (alpha, phi)
% distributed on the hemisphere with pdf (1/pi^2) sin(alpha).

% Choose hemisphere
if xi2 < 0.5
    sign_z = 1;
    xi2 = 2*xi2;
else
    sign_z = -1;
    xi2 = 2 - 2*xi2;
end

% initial guess
t0 = xi2^(2/3);

% cdf as a function of t = (1 - cos(alpha))
cdf_t = @(t) 2/pi * (acos(1-t) - sqrt(2*t-t^2) * (1 - t));

% solve for t
t = fzero(@(t) cdf_t(t) / sqrt(t) - xi2 / sqrt(t), t0);

% compute alpha and phi
alpha = acos(sign_z * (1 - t));
phi = 2*pi * xi1;

%% ---- testing code ----
N = 100000;
th = zeros(N,1);
ph = zeros(N,1);
for i = 1:N
  [th(i), ph(i)] = sin_sample(rand(), rand());
end

figure
hist(th, 100)
hold on
plot(linspace(0,pi), N * (pi/100) * 2 * sin(linspace(0,pi)).^2 / pi, 'r')
hold off

 */

MTS_IMPLEMENT_CLASS_S(DiffuseFiber, false, FiberScatteringFunction)
MTS_EXPORT_PLUGIN(DiffuseFiber, "Diffuse fiber scattering function");
MTS_NAMESPACE_END
