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
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN

class TestCore : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_atomic)
	MTS_DECLARE_TEST(test02_boundingBox)
	MTS_DECLARE_TEST(test03_boundingSphere)
	MTS_DECLARE_TEST(test04_bitmap)
	MTS_END_TESTCASE()

	void test01_atomic() {
		/* Test the atomic intrinsics in 'mitsuba/core/atomic.h' */

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

	void test02_boundingBox() {
		/* Test the generic and 3D bounding box code */
		BoundingBox3 invalid;
		assertFalse(invalid.isValid());
		assertFalse(invalid.isPoint());
		assertFalse(invalid.hasVolume());

		BoundingBox3 test(Point(1,2,3), Point(2,4,6));
		assertTrue(test == test);
		assertFalse(test != test);
		assertTrue(test != invalid);
		assertFalse(test == invalid);
		assertEquals(test.getVolume(), (Float) 6);
		assertEquals(test.getSurfaceArea(), (Float) 22);
		assertTrue(test.isValid());
		assertTrue(test.hasVolume());
		assertTrue(BoundingBox3(Point(1,2,3)).isPoint());

		assertEquals(test.getCenter(), Point(1.5f, 3.0f, 4.5f));
		assertTrue(test.contains(test.getCenter()));
		assertTrue(!test.contains(test.getCenter() + Point(10, 0, 0)));
		assertTrue(test.overlaps(test));
		assertTrue(test.contains(test));
		assertTrue(!test.overlaps(invalid));
		assertTrue(test.contains(invalid));
		assertTrue(!invalid.overlaps(test));
		assertTrue(!invalid.contains(test));
		assertEquals(test.getMajorAxis(), 2);
		assertEquals(test.getMinorAxis(), 0);

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
		assertFalse(test2.contains(test, true));
		assertFalse(test.contains(test2));
		assertTrue(test2.isValid());

		assertEquals(test.distanceTo(Point(-2, 3.5, 4.5)), (Float) 2);
		assertEquals(test.distanceTo(Point(1, 3.5, 4.5)), (Float) 0);

		assertEquals(test.distanceTo(test2), (Float) 0);
		assertEquals(test.distanceTo(BoundingBox3(Point(-2, -100, -100), 
					Point(-2, 100, 100))), (Float) 2);
		assertTrue(test2.overlaps(test2.getBoundingSphere()));

		assertFalse(test2.overlaps(BoundingSphere3(Point(-2, 3.5, 4.5), 0.5f)));

		Float nearT, farT;
		assertTrue(test2.rayIntersect(Ray(Point(0,3,4), Vector(1,0,0), 0.0f), nearT, farT));
		assertEquals(nearT, (Float) -1);
		assertEquals(farT, (Float) 2);
	
		assertFalse(test2.rayIntersect(Ray(Point(0,10,4), Vector(1,0,0), 0.0f), nearT, farT));
		assertTrue(BoundingBox3(Point(-1,-1,-1), Point(1,1,1)).getBoundingSphere() ==
				BoundingSphere3(Point(0,0,0), std::sqrt(3)));

		/* Also test the generic version of getSurfaceArea() */
		BoundingBox<Point> test3(Point(1,2,3), Point(2,4,6));
		assertEquals(test3.getSurfaceArea(), (Float) 22);

		ref<MemoryStream> mStream = new MemoryStream();
		test3.serialize(mStream);
		mStream->setPos(0);

		BoundingBox<Point> test4(mStream);
		assertTrue(test3 == test4);
	}

	void test03_boundingSphere() {
		BoundingSphere3 bsphere(Point(1,2,3), 4);
		assertTrue(bsphere.contains(Point(1, 3, 5)));
		assertTrue(bsphere.contains(Point(5, 2, 3), false));
		assertFalse(bsphere.contains(Point(5, 2, 3), true));
		assertTrue(bsphere == bsphere);
		assertFalse(bsphere != bsphere);
		bsphere.expandBy(Point(6,2,3));
		assertEquals(bsphere.getRadius(), (Float) 5);

		Float nearT = 0, farT = 0;
		assertFalse(bsphere.rayIntersect(Ray(Point(1, 6, 7), Vector(1, 0, 0), 0.0f), nearT, farT));
		assertTrue(bsphere.rayIntersect(Ray(Point(1, 2, 3), Vector(1, 0, 0), 0.0f), nearT, farT));
		assertEquals(nearT, (Float) -5);
		assertEquals(farT, (Float) 5);

		ref<MemoryStream> mStream = new MemoryStream();
		bsphere.serialize(mStream);
		mStream->setPos(0);

		BoundingSphere<Point> bsphere2(mStream);
		assertTrue(bsphere == bsphere2);
	}

	void test04_bitmap() {
		/* Test the metadata support for PNG bitmaps */
		ref<Bitmap> bitmap = new Bitmap(Vector2i(2, 2), 8);
		uint8_t *data = bitmap->getData();
		for (int i=0; i<4; ++i)
			data[i] = i;
		bitmap->setGamma(1.75f);
		bitmap->setTitle("Test title");
		bitmap->setAuthor("Mitsuba renderer");
		bitmap->setComment("A comment");

		ref<MemoryStream> fs = new MemoryStream();
		bitmap->save(Bitmap::EPNG, fs);
		fs->setPos(0);

		ref<Bitmap> bitmap2 = new Bitmap(Bitmap::EPNG, fs);
		uint8_t *data2 = bitmap2->getData();
		for (int i=0; i<4; ++i)
			assertEquals(data2[i], i);
		assertEquals(bitmap2->getGamma(), 1.75f);
		assertEquals(bitmap2->getTitle(), "Test title");
		assertEquals(bitmap2->getAuthor(), "Mitsuba renderer");
		assertEquals(bitmap2->getComment(), "A comment");
	}
};

MTS_EXPORT_TESTCASE(TestCore, "Testcase for core routines")
MTS_NAMESPACE_END
