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

#if !defined(__TESTCASE_H)
#define __TESTCASE_H

#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

#if defined(MTS_TESTCASE)

/**
 * When a testcase is being compiled, define the following preprocessor macros for convenience
 */
#define assertEquals(actual, expected) assertEqualsImpl(actual, expected, 0, __FILE__, __LINE__)
#define assertEqualsEpsilon(actual, expected, epsilon) assertEqualsImpl(actual, expected, epsilon, __FILE__, __LINE__)
#define assertTrue(expr) assertTrueImpl(expr, #expr, __FILE__, __LINE__)
#define assertFalse(expr) assertFalseImpl(expr, #expr, __FILE__, __LINE__)
#endif

/** \brief Base class of all testcases. Implementations of this
 * interface can be executed using the 'mtsutil' command. The execution
 * order is as follows: after initializaiton using init(), any tests
 * declared using the MTS_DECLARE_TEST() macro are executed. Finally,
 * the shutdown() method is called. See the files in 'mitsuba/src/tests'
 * for examples.
 */
class MTS_EXPORT_RENDER TestCase : public Utility {
public:
	/**
	 * Perform any required initializations. The default
	 * implementation simply returns
	 */
	virtual void init();

	/**
	 * Execute any required shutdown code. The default
	 * implementation simply returns
	 */
	virtual void shutdown();

	/// Return the number of executed testcases
	inline int getExecuted() const { return m_executed; }
	
	/// Return the number of successfully executed testcases
	inline int getSucceeded() const { return m_succeeded; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~TestCase() { }

	/// Asserts that the two integer values are equal
	void assertEqualsImpl(int actual, int expected, Float, const char *file, int line);

	/// Asserts that the two string values are equal
	void assertEqualsImpl(const std::string &actual, const std::string &expected, Float, const char *file, int line);
	
	/// Asserts that the two single precision floating point values are equal
	void assertEqualsImpl(double actual, double expected, Float epsilon, const char *file, int line);
	
	/// Asserts that the two double precision floating point values are equal
	void assertEqualsImpl(float actual, float expected, Float epsilon, const char *file, int line);

	/// Asserts that the two fixed-size or dynamic matrices are equal
	template<typename T, int M1, int N1, int M2, int N2> 
			void assertEqualsImpl(const Eigen::Matrix<T, M1, N1> &actual,
				const Eigen::Matrix<T, M2, N2> &expected, Float epsilon, const char *file, int line) {
		bool match = true;
		if (expected.rows() != actual.rows() || expected.cols() != actual.cols())
			Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
				"mismatched matrix dimensions -- expected %ix%i, got %ix%i.", 
				(int) expected.rows(), (int) expected.cols(), (int) actual.rows(), (int) actual.cols());
		for (int i=0; i<expected.rows(); ++i)
			for (int j=0; j<expected.cols(); ++j)
				if (std::abs(expected(i, j)-actual(i, j)) > epsilon)
					match = false;
		if (!match)
			Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
				"expected matrix %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
	}

	/// Asserts that a condition is true
	void assertTrueImpl(bool condition, const char *expr, const char *file, int line);

	/// Asserts that a condition is false
	void assertFalseImpl(bool condition, const char *expr, const char *file, int line);
protected:
	int m_executed, m_succeeded;
};

MTS_NAMESPACE_END

#define EXECUTE_GUARDED(name) \
	try { \
		Log(EInfo, "Executing test \"%s\" ..", #name); \
		m_executed++;\
		name();\
		m_succeeded++;\
	} catch (std::exception &e) {\
		Log(EInfo, "Testcase failed with error: %s", e.what());\
	}

#define MTS_BEGIN_TESTCASE() \
	MTS_DECLARE_CLASS() \
	int run(int argc, char **argv) {\
		init(); \
		Log(EInfo, "Executing testcase \"%s\" ..", getClass()->getName().c_str()); \
		m_executed = m_succeeded = 0;

#define MTS_DECLARE_TEST(name) \
		EXECUTE_GUARDED(name)

#define MTS_END_TESTCASE()\
		shutdown();\
		return 0;\
	}

#define MTS_EXPORT_TESTCASE(name, descr) \
	MTS_IMPLEMENT_CLASS(name, false, TestCase) \
	extern "C" { \
		void MTS_EXPORT *CreateUtility() { \
			return new name(); \
		} \
		const char MTS_EXPORT *GetDescription() { \
			return descr; \
		} \
	}

#endif /* __TESTCASE_H */
