/*
 * bitslice.h
 *
 *  Created on: 5 Feb 2016
 *      Author: lennonjo
 */

#ifndef BITSLICE_H_
#define BITSLICE_H_

class bitslice {
public:
	bitslice();
	virtual ~bitslice();

	unsigned long long __rdtsc (void);

	template <typename T, typename P>
	inline void bitslice_convert (P * numbers, T * bitslices, const int nSize);

	template <typename T, typename P>
	inline void bitslice_convert_back (T * bitslices, P * numbers, const int bSize);

	template <typename T>
	inline void bitslice_shift_left(T a[], const int size, const int numShifts);

	template <typename T>
	inline void bitslice_shift_right(T a[], const int size, const int numShifts);

	template <typename T>
	inline void bitslice_arithmetic_shift_right(T a[], const int size, const int numShifts);

	template <typename T>
	inline void bitslice_add(T * __restrict__ a, T * __restrict__ b, T c[], const int size);

	template <typename T>
	inline void bitslice_sub(T * __restrict__ a, T * __restrict__ b, T c[], const int size);

	template <typename T>
	inline void bitslice_mul(T * __restrict__  a, T * __restrict__  b, T * __restrict__ c, const int size);

	template <typename T>
	inline void bitslice_signed_mul(T * __restrict__ a, T * __restrict__ b, T * __restrict__ c, const int size);

	template <typename T>
	inline void bitslice_fp_signed_mul(T * __restrict__ a, T * __restrict__ b, T * __restrict__ c, const int size, const int places);

	template <typename T>
	inline void bitslice_horizontal_add(T * __restrict__ a, T * __restrict__ result, const int size);

	template <typename T>
	inline void bitslice_blas_scale(T * __restrict__ vec, T * __restrict__ scale, T * __restrict__ result, const int size);

	template <typename T>
	inline void bitslice_blas_dot(T * __restrict__ a, T * __restrict__ b, T * __restrict__ c, const int size);

	template <typename T>
	inline void insert_random_values (T * a, const int limit, const int size);

	template <typename T>
	inline void print_bitslice_arrays(T a[], T b[], T c[], const int size);

	template <typename T>
	inline void print_bitslice_arrays_demo(T a[], T b[], T c[], const int size);

	template <typename T>
	inline void print_bitslice_array(T a[], const int size);

	template <typename T>
	inline void print_bitslice_fp_array(T a[], const int size, const int decimal_place);

	template <typename T>
	inline void print_clock_cycles(T * cycles, const int length);
};

#endif /* BITSLICE_H_ */
