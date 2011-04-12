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

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/quad.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

class TestQuadrature : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_quad)
	MTS_DECLARE_TEST(test02_nD_01)
	MTS_DECLARE_TEST(test03_nD_02)
	MTS_END_TESTCASE()

	Float testF(Float t) const {
		return std::sin(t);
	}

	void testF2(const Float *in, Float *out) const {
		*out = std::sin(*in);
	}
	
	inline Float gauss3(Vector x, Float stddev) const {
		return std::exp(-0.5f * dot(x, x)/stddev)/(std::pow(2*M_PI * stddev, (Float) 3 / (Float) 2));
	}

	void testF3(size_t nPoints, const Float *in, Float *out) const {
		for (size_t i=0; i<nPoints; ++i) {
			Vector v(in[0], in[1], in[2]);
			Float weight = 
				(1 + v.x*v.x) / std::pow(1-v.x*v.x, 2) *
				(1 + v.y*v.y) / std::pow(1-v.y*v.y, 2) *
				(1 + v.z*v.z) / std::pow(1-v.z*v.z, 2);
			v = Vector(v.x / (1-v.x*v.x), v.y / (1-v.y*v.y),
					   v.z / (1-v.z*v.z));
			out[0*nPoints + i] = weight * gauss3(v, 0.1);
			out[1*nPoints + i] = weight * gauss3(v, 1);
			in += 3;
		}
	}

	void test01_quad() {
		GaussLobattoIntegrator quad(1024, 0, 1e-5f);
		size_t evals;
		Float result = quad.integrate(boost::bind(
			&TestQuadrature::testF, this, _1), 0, 10, evals);
		Float ref = 2 * std::pow(std::sin(5), 2);
		Log(EInfo, "test01_quad(): used " SIZE_T_FMT " function evaluations", evals);
		assertEqualsEpsilon(result, ref, 1e-5f);
	}

	void test02_nD_01() {
		NDIntegrator quad(1, 1, 1024, 0, 1e-5f);
		Float min = 0, max = 10, result, err;
		size_t evals;
		assertEquals(quad.integrate(boost::bind(
			&TestQuadrature::testF2, this, _1, _2), &min, &max, &result, &err, evals),
			NDIntegrator::ESuccess);
		Float ref = 2 * std::pow(std::sin(5), 2);
		Log(EInfo, "test02_nD_01(): used " SIZE_T_FMT " function evaluations, "
				"error=%f", evals, err);
		assertEqualsEpsilon(result, ref, 1e-5f);
	}

	void test03_nD_02() {
		NDIntegrator quad(2, 3, 1000000, 0, 1e-5f);
		size_t evals;
		Float min[3] = { -1, -1, -1 } , max[3] = { 1, 1, 1 }, result[2], err[2];
		assertEquals(quad.integrateVectorized(boost::bind(
			&TestQuadrature::testF3, this, _1, _2, _3), min, max, result, err, evals),
			NDIntegrator::ESuccess);
		Log(EInfo, "test02_nD_02(): used " SIZE_T_FMT " function evaluations, "
				"error=[%f, %f]", evals, err[0], err[1]);
		assertEqualsEpsilon(result[0], 1, 1e-5f);
		assertEqualsEpsilon(result[1], 1, 1e-5f);
	}
};

MTS_EXPORT_TESTCASE(TestQuadrature, "Testcase for quadrature routines")
MTS_NAMESPACE_END
