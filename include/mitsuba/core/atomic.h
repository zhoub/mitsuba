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

#if !defined(__ATOMIC_H)
#define __ATOMIC_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

#if defined(__OSX__)
#include <libkern/OSAtomic.h>
#endif

#if defined(WIN32)
#include <intrin.h>
#endif

/**
 * \brief Implements a range of atomic operations on builtin data types
 */
class Atomic {
public:
	/**
	 * \brief Atomically attempt to exchange a pointer with another value
	 * 
	 * \param v Pointer to the pointer in question
	 * \param oldValue Last known value of the destination \a v
	 * \param newValue Replacement value for the destination \a v
	 * \tparam T Base type of the pointer
	 * \return \c true if \c *v was equal to \c oldValue and the exchange
	 *         was successful.
	 */
	template <typename T> inline static 
		bool compareAndExchangePtr(T **v, T *newValue, T *oldValue) {
	#if defined(WIN32)
	    return InterlockedCompareExchangePointer(
			reinterpret_cast<volatile PVOID *>(v), newValue, oldValue) == oldValue;
	#else
	  #if !defined(__clang__)
		return __sync_bool_compare_and_swap(v, oldValue, newValue);
	  #else
	  #if __SIZEOF_POINTER__ == 8 
		return __sync_bool_compare_and_swap(
			reinterpret_cast<long long volatile *>(v), reinterpret_cast<long long>(oldValue), 
			reinterpret_cast<long long>(newValue));
	  #else
		return __sync_bool_compare_and_swap(
			reinterpret_cast<long volatile *>(v), reinterpret_cast<long>(oldValue), 
			reinterpret_cast<long>(newValue));
	  #endif
	  #endif
	#endif
	}
	
	/**
	 * \brief Atomically attempt to exchange a 32-bit integer with another value
	 * 
	 * \param v Pointer to the memory region in question
	 * \param oldValue Last known value of the destination \a v
	 * \param newValue Replacement value for the destination \a v
	 * \return \c true if \c *v was equal to \c oldValue and the exchange
	 *         was successful.
	 */
	inline static bool compareAndExchange(volatile int32_t *v, int32_t newValue, int32_t oldValue) {
	#if defined(WIN32)
	    return InterlockedCompareExchange(
			reinterpret_cast<volatile LONG *>(v), newValue, oldValue) == oldValue;
	#else
		return __sync_bool_compare_and_swap(v, oldValue, newValue);
	#endif
	}
	
	/**
	 * \brief Atomically attempt to exchange a 64-bit integer with another value
	 * 
	 * \param v Pointer to the memory region in question
	 * \param oldValue Last known value of the destination \a v
	 * \param newValue Replacement value for the destination \a v
	 * \return \c true if \c *v was equal to \c oldValue and the exchange
	 *         was successful.
	 */
	
	inline static bool compareAndExchange(volatile int64_t *v, int64_t newValue, int64_t oldValue) {
	#if defined(WIN32)
	    return _InterlockedCompareExchange64(
			reinterpret_cast<volatile LONGLONG *>(v), newValue, oldValue) == oldValue;
	#else
		return __sync_bool_compare_and_swap(v, oldValue, newValue);
	#endif
	}
	
	/**
	 * \brief Atomically add \a delta to the 32-bit integer destination \a dst
	 *
	 * \return The final value written to \a dst
	 */
	inline static int32_t add(volatile int32_t *dst, int32_t delta) {
		#if defined(WIN32)
			return InterlockedExchangeAdd(dst, delta) + delta;
		#else
			return __sync_add_and_fetch(dst, delta);
		#endif
	}
	
	/**
	 * \brief Atomically add \a delta to the 64-bit integer destination \a dst
	 *
	 * \return The final value written to \a dst
	 */
	inline static int64_t add(volatile int64_t *dst, int64_t delta) {
		#if defined(WIN64)
			return _InterlockedExchangeAdd64(dst, delta) + delta;
		#elif defined(WIN32)
			SLog(EError, "add() cannot handle 64-bit integers on WIN32");
			return 0;
		#else
			return __sync_add_and_fetch(dst, delta);
		#endif
	}
	
	/**
	 * \brief Atomically add \a delta to the floating point destination \a dst
	 *
	 * \return The final value written to \a dst
	 */
	inline static float add(volatile float *dst, float delta) {
		/* Atomic FP addition from PBRT */
		union bits { float f; int32_t i; };
		bits oldVal, newVal;
		do {
			// From PBRT:
	        // On IA32/x64, adding a PAUSE instruction in compare/exchange loops
	        // is recommended to improve performance.  (And it does!)
			#if (defined(__i386__) || defined(__amd64__))
	        	__asm__ __volatile__ ("pause\n");
			#endif
			oldVal.f = *dst;
			newVal.f = oldVal.f + delta;
		} while (!compareAndExchange((volatile int32_t *) dst, newVal.i, oldVal.i));
		return newVal.f;
	}
	
	/**
	 * \brief Atomically add \a delta to the floating point destination \a dst
	 *
	 * \return The final value written to \a dst
	 */
	inline static double add(volatile double *dst, double delta) {
		/* Atomic FP addition from PBRT */
		union bits { double f; int64_t i; };
		bits oldVal, newVal;
		do {
			// From PBRT:
	        // On IA64/x64, adding a PAUSE instruction in compare/exchange loops
	        // is recommended to improve performance.  (And it does!)
			#if (defined(__i386__) || defined(__amd64__))
		        __asm__ __volatile__ ("pause\n");
			#endif
			oldVal.f = *dst;
			newVal.f = oldVal.f + delta;
		} while (!compareAndExchange((volatile int64_t *) dst, newVal.i, oldVal.i));
		return newVal.f;
	}
};

MTS_NAMESPACE_END
	
#endif /* __ATOMIC_H */
