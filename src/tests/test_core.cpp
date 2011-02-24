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

#include <mitsuba/core/testcase.h>
#include <mitsuba/core/bbox.h>
#include <mitsuba/core/atomic.h>

MTS_NAMESPACE_BEGIN

class TestCore : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_boundingBox)
	MTS_DECLARE_TEST(test02_atomic)
	MTS_END_TESTCASE()

	void test01_boundingBox() {
		BoundingBox3 invalid;
		assertFalse(invalid.isValid());

		BoundingBox3 test(Point(1,2,3), Point(2,4,6));
		assertTrue(test == test);
		assertFalse(test == invalid);
		assertEquals(test.getVolume(), 6);
		assertEquals(test.getSurfaceArea(), 22);
		assertTrue(test.isValid());

		assertEquals(test.getCenter(), Point(1.5f, 3.0f, 4.5f));
		assertTrue(test.contains(test.getCenter()));
		assertTrue(!test.contains(test.getCenter() + Point(10, 0, 0)));
		assertTrue(test.overlaps(test));
		assertTrue(test.contains(test));
		assertTrue(!test.overlaps(invalid));
		assertTrue(!test.contains(invalid));
		assertTrue(!invalid.overlaps(test));
		assertTrue(!invalid.contains(test));
		assertEquals(test.getLongestDimension(), 2);
		assertEquals(test.getShortestDimension(), 0);

		test.expandBy(Point(0,5,6));
		assertTrue(test == BoundingBox3(Point(0, 2, 3), Point(2,5,6)));

		BoundingBox3 test2(Point(-1, 2, 3), Point(1, 5, 6));
		assertTrue(test2.overlaps(test));
		assertTrue(test.overlaps(test2));

		assertTrue(test2.overlaps(test));
		assertTrue(test.overlaps(test2));
		assertFalse(test2.contains(test));
		assertFalse(test.contains(test2));

		test2.expandBy(test);
		assertTrue(test2.contains(test));
		assertFalse(test.contains(test2));
		assertTrue(test2.isValid());

		assertEquals(test.distanceTo(Point(-2, 3.5, 4.5)), 1);
		assertEquals(test.distanceTo(Point(1, 3.5, 4.5)), 0);

		assertEquals(test.distanceTo(test2), 0);
		assertEquals(test.distanceTo(BoundingBox3(Point(-2, -100, -100), 
					Point(-2, 100, 100))), 1);
		assertTrue(test2.overlaps(test2.getBoundingSphere()));

		assertFalse(test2.overlaps(BoundingSphere(Point(-2, 3.5, 4.5), 0.5f)));
	
		Float nearT, farT;
		assertTrue(test2.rayIntersect(Ray(Point(0,0,0), Vector(1,0,0), 0.0f), nearT, farT));
		assertEquals(nearT, -1);
		assertEquals(farT, 2);
	
		assertFalse(test2.rayIntersect(Ray(Point(0,10,0), Vector(1,0,0), 0.0f), nearT, farT));

		/* Also test the generic version of getSurfaceArea() */
		BoundingBox<Point> test3(Point(1,2,3), Point(2,4,6));
		assertEquals(test3.getSurfaceArea(), 22);
	}
	
	void test02_atomic() {
		int A = 0, B = 1, *C = &A;
		int64_t A64 = 0;

		assertTrue(Atomic::compareAndExchangePtr(&C, &B, &A));
		assertTrue(C == &B);

		assertTrue(Atomic::compareAndExchange(&A, 5, 0));
		assertEquals(A, 5);
		
		assertTrue(Atomic::compareAndExchange(&A64, 5, 0));
		assertEquals(A64, 5);

		assertEquals(Atomic::add(&A, 1), 6);
		assertEquals(A, 6);
		assertEquals(Atomic::add(&A64, 1), 6);
		assertEquals(A64, 6);

		float Af = 5.0f;
		double Ad = 5.0;
		assertEquals(Atomic::add(&Af, 1), 6.0f);
		assertEquals(Af, 6.0f);
		assertEquals(Atomic::add(&Ad, 1), 6.0);
		assertEquals(Ad, 6.0);
	}
};

MTS_EXPORT_TESTCASE(TestCore, "Testcase for core routines")
MTS_NAMESPACE_END
