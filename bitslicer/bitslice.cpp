//#include "bitslice.h"

//bitslice::bitslice() {
	// TODO Auto-generated constructor stub

//}

//bitslice::~bitslice() {
	// TODO Auto-generated destructor stub
//}

#include <iostream>
#include <inttypes.h>
#include <stdio.h>
#include <bitset>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
using namespace std;

unsigned long long __rdtsc (void) {
  unsigned cycles_low, cycles_high;
  asm volatile ("CPUID\n\t"
		"RDTSC\n\t"
		"mov %%edx, %0\n\t"
		"mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::
		"%rax", "%rbx", "%rcx", "%rdx");
  return ((unsigned long long)cycles_high << 32) | cycles_low;
}

/* Convert an array of integers into a bitslice array */
template <typename T, typename P>
inline void bitslice_convert (P * numbers, T * bitslices, const int nSize)
{
    if (nSize != sizeof(T)*8)
    {
        std::cout << "CONVERSION ERROR: the amount of integers in 'numbers' must match the width of the bitslice container.\nDetails;\n";
        std::cout << "numbers length: " << nSize << std::endl;
        std::cout << "bitslice width: " << sizeof(T)*8 << std::endl;
    }

    T mask = 1;
    T element = 0;
    int n = 0;

    for (int i = 0; i < sizeof(P)*8; i++)
    {
        for (int j = 0; j < nSize; j++)
        {
            element |= ( (numbers[j] & mask) >> i ) << n++;
        }
        bitslices[i] = element;

        element = 0;
        mask <<= 1;
        n = 0;
    }
}

/* Convert a bitslice array back into a conventional array of integers */
template <typename T, typename P>
inline void bitslice_convert_back (T * bitslices, P * numbers, const int bSize)
{
    if (bSize != sizeof(P)*8)
    {
        std::cout << "CONVERSION ERROR: the length of the bitslice array must match the width of the 'numbers' type.\nDetails;\n";
        std::cout << "bitslice length: " << bSize << std::endl;
        std::cout << "numbers width: " << sizeof(P)*8 << std::endl;
    }

    T mask = 1;
    P element = 0;

    for (int i = 0; i < sizeof(T)*8; i++)
    {
        for (int j = 0; j < sizeof(P)*8; j++)
        {
            element |= ( (bitslices[j] & mask) >> i ) << j;
        }
        numbers[i] = element;
        element = 0;
        mask <<= 1;
    }
}

/*
    Left-shift a bitslice structure by 'numShifts'.
*/
template <typename T>
inline void bitslice_shift_left(T a[], const int size, const int numShifts)
{
    if (numShifts > 0)
    {
        int count = 0;
        while (count < numShifts)
        {
            for (int i = size-1; i >= 1; i--)
            {
                a[i] = a[i-1];
            }
            a[0] = 0;
            count++;
        }
    }

    /* Alternative (but slower) implementation */
    //memmove(a+1, a, sizeof(T)*(size-1));
    //a[0] = 0;
}

/*
    Right-shift a bitslice structure by 'numShifts'.
*/
template <typename T>
inline void bitslice_shift_right(T a[], const int size, const int numShifts)
{
    if (numShifts > 0)
    {
        int count = 0;
        while (count < numShifts)
        {
            for (int i = 0; i < size-1; i++)
            {
                a[i] = a[i+1];
            }
            a[size-1] = 0;
            count++;
        }
    }

    /* Alternative (but slower) implementation */
    //memmove(a, a+1, sizeof(T)*(size-1));
    //a[size-1] = 0;
}

/*
    Arithmetic right-shift by 'numShifts'.
*/
template <typename T>
inline void bitslice_arithmetic_shift_right(T a[], const int size, const int numShifts)
{
    if (numShifts > 0)
    {
        int count = 0;
        while (count < numShifts)
        {
            for (int i = 0; i < size-1; i++){
                a[i] = a[i+1];
            }
            a[size-1] = a[size-2];
            count++;
        }
    }

    /* Alternative (but slower) implementation */
    //T tmp = a[size-2];
    //memmove(a, a+1, sizeof(T)*(size-1));
    //a[size-1] = tmp;
}

/*
    This function adds all bitsliced values contained in B
    to all bitsliced values contained in A, storing all bitsliced
    results in C. Functionality mapped from gate-level logic.

    Gate Logic: http://upload.wikimedia.org/wikipedia/commons/6/69/Full-adder_logic_diagram.svg

    @params size - the number of bits in each bitsliced value.
*/
template <typename T>
inline void bitslice_add(T * __restrict__ a, T * __restrict__ b, T c[], const int size)
{
    /*T carry = 0;
    T xxor = a[0] ^ b[0];
    T aand = a[0] & b[0];

    c[0] = carry ^ xxor;
    carry = aand;

    for (int i = 1; i < size; i++)
    {
        xxor = a[i] ^ b[i];
        aand = a[i] & b[i];
        c[i] = carry ^ xxor;
        carry = (carry & xxor) | aand;
    }*/

    c[0] = a[0] ^ b[0];
    T carry = a[0] & b[0];

    for (int i = 1; i < size; i++)
    {
        c[i] = carry ^ (a[i] ^ b[i]);
        carry = (carry & (a[i] ^ b[i])) | (a[i] & b[i]);
    }
}

/*
    This function subtracts all bitsliced values contained in B
    from all bitsliced values contained in A, storing all bitsliced
    results in C. Functionality mapped from gate-level logic.

    Gate Logic: http://ustudy.in/sites/default/files/fullsub.gif

    @params size - the number of bits in each bitsliced value.
*/
template <typename T>
inline void bitslice_sub(T * __restrict__ a, T * __restrict__ b, T c[], const int size)
{
    T borrow = 0;
    T xxor = a[0] ^ b[0];
    c[0] = borrow ^ xxor;
    T aand1 = ~a[0] & b[0];
    T aand2 = borrow & ~xxor;
    borrow = aand1 | aand2;

    for (int i = 1; i < size; i++)
    {
        xxor = a[i] ^ b[i];
        c[i] = borrow ^ xxor;
        aand1 = ~a[i] & b[i];
        aand2 = borrow & ~xxor;
        borrow = aand1 | aand2;
    }
}

/*
    Multiply all corresponding bitslice values from A and B. Store result in C.
*/
template <typename T>
inline void bitslice_mul(T * __restrict__  a, T * __restrict__  b, T * __restrict__ c, const int size)
{
    for (int i = 0; i < size; i++)
    {
        c[i] = a[i] & b[0];
    }

    for (int i = 1; i < size; i++)
    {
        T carry = 0;
        T xxor, aand;

        for (int j = i; j < size; j++)
        {
            T current = a[j-i] & b[i];
            xxor = current ^ c[j];
            aand = current & c[j];
            c[j] = carry ^ xxor;
            carry = (carry & xxor) | aand;
        }
    }
}

/*
    A bitsliced implementation of the Baugh-Wooley signed multiplication algorithm.
    'c' must be twice the size of 'a' and 'b' in order to fully store the product of 'a' and 'b'.
    'size' still represents the sizes of 'a' and 'b' respectively.
*/
template <typename T>
inline void bitslice_signed_mul(T * __restrict__ a, T * __restrict__ b, T * __restrict__ c, const int size)
{
    /* first partial */
    for (int i = 0; i < size-1; i++)
    {
        c[i] = a[i] & b[0];
    }

    /* negate MSB of first partial */
    c[size-1] = ~(a[size-1] & b[0]);

    /* append the value 1 to the MSB of the first partial,
    for every bitslice (hence the need for all 1's) */
    c[size] = -1;

    /* intermediate partials */
    T * partial = (T *) calloc(size*2, sizeof(T));

    for (int i = 1; i < size-1; i++)
    {
        /* (old) zero padding */
        //for (int j = 0; j < i; j++) {
        //    partial[j] = 0;
        //}
        /* (new) zero padding */
        partial[i-1] = 0;

        for (int j = i; j < (i+size-1); j++)
        {
            partial[j] = a[j-i] & b[i];
        }

        partial[i+size-1] = ~(a[size-1] & b[i]);

        bitslice_add(c, partial, c, size*2);
    }

    /* final partial */

    /* (old) zero padding */
    //for (int j = 0; j < size-1; j++) {
    //    partial[j] = 0;
    //}

    /* (new) zero padding */
    partial[size-2] = 0;

    for (int j = size-1; j < size*2-2; j++)
    {
        partial[j] = ~(a[j-size-1] & b[size-1]);
    }

    partial[size*2-2] = a[size-1] & b[size-1];

    /* append 1 to the MSB of the last partial,
    for every bitslice (hence the need for all 1's) */
    partial[size*2-1] = -1;

    bitslice_add(c, partial, c, size*2);

    free(partial);
}

/*
    Fixed-point signed multiplication (signed multiplication followed by an arithmetic right shift).
    Both operands must have the same number of decimal places. The initial product has twice the
    number of decimal places, hence the need to shift the product by the original number of places.
*/
template <typename T>
inline void bitslice_fp_signed_mul(T * __restrict__ a, T * __restrict__ b, T * __restrict__ c, const int size, const int places)
{
    bitslice_signed_mul(a,b,c,size);
    bitslice_arithmetic_shift_right(c,size*2,places);
}

/*
    Sum all bitslice values contained in 'a' and
    store result in leftmost bitslice position of 'result'.
*/
template <typename T>
inline void bitslice_horizontal_add(T * __restrict__ a, T * __restrict__ result, const int size)
{
    T * tmp = (T *) malloc (sizeof(T) * size);

    for (int i = 0; i < size; i++)
    {
        tmp[i] = a[i] << 1;
    }

    bitslice_add(a,tmp,result,size);

    int count = 0;

    while (count < sizeof(T)*8 - 1)
    {
        for (int i = 0; i < size; i++)
        {
            tmp[i] <<= 1;
        }
        bitslice_add(result,tmp,result,size);
        count++;
    }
    free(tmp);
}

/*
    Vector Scale
*/
template <typename T>
inline void bitslice_blas_scale(T * __restrict__ vec, T * __restrict__ scale, T * __restrict__ result, const int size)
{
    bitslice_signed_mul(vec,scale,result,size);
}

/*
    Dot product
*/
template <typename T>
inline void bitslice_blas_dot(T * __restrict__ a, T * __restrict__ b, T * __restrict__ c, const int size)
{
    T * tmp = (T *) calloc (size*2,sizeof(T));
    bitslice_signed_mul(a,b,tmp,size);
    bitslice_horizontal_add(tmp,c,size*2);
    free(tmp);
}

/*
    Insert random values into array.
*/
template <typename T>
inline void insert_random_values (T * a, const int limit, const int size)
{
    for (int i = 0; i < size; i++)
    {
        a[i] = rand() % limit;
    }
}

/* Print bitslice structures */
template <typename T>
inline void print_bitslice_arrays(T a[], T b[], T c[], const int size)
{
    //for (int i = size-1; i >= 0; i--)
    for (int i = 0; i < size; i++)
    {
        std::cout << "bit " << i << ":\t"
            << (std::bitset<sizeof(T)*8>) a[i] << "\t"
            << (std::bitset<sizeof(T)*8>) b[i] << "\t"
            << (std::bitset<sizeof(T)*8>) c[i] <<
        std::endl;
    }
}

template <typename T>
inline void print_bitslice_arrays_demo(T a[], T b[], T c[], const int size)
{
    for (int i = 0; i < size; i++)
    //for (int i = size-1; i >= 0; i--)
    {
        std::cout << "bit " << i << ":\t"
            << (std::bitset<sizeof(T)*8>) a[i] << "\t"
            << (std::bitset<sizeof(T)*8>) b[i] << "\t"
            << std::endl;
    }

    std::cout << "\nResult;\n";
    for (int i = 0; i < size; i++)
    //for (int i = size-1; i >= 0; i--)
    {
        std::cout << "bit " << i << "\t" << (std::bitset<sizeof(T)*8>) c[i] << std::endl;
    }
}

/* Print single bitslice structure */
template <typename T>
inline void print_bitslice_array(T a[], const int size)
{
    //for (int i = size-1; i >= 0; i--)
    for (int i = 0; i < size; i++)
    {
        std::cout << "bit " << i << ":\t" << (std::bitset<sizeof(T)*8>) a[i] << std::endl;
    }
}

/* Print fixed-point bitslice structure */
template <typename T>
inline void print_bitslice_fp_array(T a[], const int size, const int decimal_place)
{
    std::cout << "Fixed-Point Output (" << decimal_place << " decimal place/s);" << std::endl;
    //for (int i = size-1; i >= 0; i--)
    for (int i = 0; i < size; i++)
    {
        if (i == decimal_place)
        {
            switch(sizeof(T)*8)
            {
                case 64:
                    std::cout << "dec  :\t................................................................" << std::endl;
                    break;
                case 32:
                    std::cout << "dec  :\t................................" << std::endl;
                    break;
                case 16:
                    std::cout << "dec  :\t................" << std::endl;
                    break;
                case 8:
                    std::cout << "dec  :\t........" << std::endl;
                    break;
            }

            //std::cout << "dec  :\t................................................................" << std::endl;
            std::cout << "bit " << i << ":\t" << (std::bitset<sizeof(T)*8>) a[i] << std::endl;
        }
        else
        {
            std::cout << "bit " << i << ":\t" << (std::bitset<sizeof(T)*8>) a[i] << std::endl;
        }
    }
}

/* Source: http://xoax.net/cpp/ref/cpp_examples/incl/mean_med_mod_array/ */
double getMedian(uint64_t daArray[], int iSize)
{
    // Allocate an array of the same size and sort it.
    uint64_t * dpSorted;
    uint64_t dTemp;

    dpSorted = new uint64_t[iSize];

    for (int i = 0; i < iSize; ++i) {
        dpSorted[i] = daArray[i];
    }

    for (int i = iSize - 1; i > 0; --i)
    {
        for (int j = 0; j < i; ++j)
        {
            if (dpSorted[j] > dpSorted[j+1])
            {
                dTemp = dpSorted[j];
                dpSorted[j] = dpSorted[j+1];
                dpSorted[j+1] = dTemp;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    double dMedian = 0.0;
    if ((iSize % 2) == 0)
    {
        dMedian = (dpSorted[iSize/2] + dpSorted[(iSize/2) - 1])/2.0;
    }
    else
    {
        dMedian = dpSorted[iSize/2];
    }
    delete [] dpSorted;
    return dMedian;
}

template <typename T>
inline void print_clock_cycles(T * cycles, const int length)
{
    std::cout << "Print clock cycles\n";
    for (int i = 0; i < length; i++){
        std::cout << cycles[i] << ",";
    }
}

typedef uint64_t base_type;

//typedef uint64_t std_type;
//typedef uint32_t std_type;
//typedef uint16_t std_type;
typedef uint8_t std_type;

int main()
{
    const int PASSES = 100;

    srand (time(NULL));

    base_type start, stop;

    base_type add_timings[PASSES];
    base_type sub_timings[PASSES];
    base_type mul_timings[PASSES];
    base_type signed_mul_timings[PASSES];
    base_type fixed_point_timings[PASSES];
    base_type scale_timings[PASSES];
    base_type horizontal_add_timings[PASSES];
    base_type dot_timings[PASSES];
    base_type normal_add_times[PASSES];
    base_type normal_sub_times[PASSES];
    base_type normal_mul_times[PASSES];
    base_type normal_hAdd_times[PASSES];
    base_type normal_vector_scale_timings[PASSES];
    base_type normal_dot_product_timings[PASSES];
    base_type convert_timings[PASSES];
    base_type convert_back_timings[PASSES];

    /* The number of bits contained in each bitslice value */
    const int size = 8;

    std::cout << "BITSLICE VECTOR COMPUTATION\n\n";
    std::cout << "64 x " << size << "-bit bitslice values.\n\n";
    std::cout << "Note: each bitslice is to be read vertically.\n";

    /* Enter values that will translate to appropriate bitslice values */
    /* 8-bit */
    base_type a[] = {-1,-1,1,0,0,0,0,0};
    base_type b[] = {3,-1,1,0,0,0,0,0};

    /* 16-bit */
    //base_type a[] = {-1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //base_type b[] = {3,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};

    /* 32-bit */
    //base_type a[] = {-1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //base_type b[] = {3,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    /* 64-bit */
    //base_type a[] = {-1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //base_type b[] = {3,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    /* Used in the bitslice scale function. Scales A by the constant 3 */
    base_type scale[] = {-1,-1,0,0,0,0,0,0};
    //base_type scale[] = {-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //base_type scale[] = {-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //base_type scale[] = {-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //base_type scale[] = {-1,-1,0,0,0,0,0,0,0};

    /* Normal values (matching those from the bitsliced structure) */
    std_type alpha[] = {7,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
    std_type beta[] = {7,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};

    /* Storage for bitslice results */
    base_type * sums = (base_type *) calloc (size, sizeof(base_type));
    base_type * diffs = (base_type *) calloc (size, sizeof(base_type));
    base_type * products = (base_type *) calloc (size, sizeof(base_type));
    base_type * productsDouble = (base_type *) calloc (size*2, sizeof(base_type));
    base_type * signed_result = (base_type *) calloc(size*2,sizeof(base_type));
    base_type * fixed_point_result = (base_type *) calloc (size*2, sizeof(base_type));
    base_type * horizontal_add_result = (base_type *) calloc (size, sizeof(base_type));
    base_type * dot_result = (base_type *) calloc (size*2, sizeof(base_type));
    base_type * scale_result = (base_type *) calloc(size*2,sizeof(base_type));

    /* Storage for conventional results */
    std_type * normal_add_results = (std_type *) calloc (sizeof(base_type)*8,sizeof(std_type));
    std_type * normal_sub_results = (std_type *) calloc (sizeof(base_type)*8,sizeof(std_type));
    std_type * normal_mul_results = (std_type *) calloc (sizeof(base_type)*8,sizeof(std_type));
    std_type * normal_vector_scale_result = (std_type *) calloc (sizeof(base_type)*8,sizeof(std_type));
    int normal_hAdd_result = 0;
    int normal_dot_product_result = 0;

    /* containers for clock cycle timings */
    int uAdd, uSub, uMul, sMul; // bitslice add/sub/mul/sMul
    int fixp;                   // bitslice fixed-point mul
    int hAdd;                   // bitslice horizontal adder
    int bDot;                   // bitslice dot product
    int scal;                   // bitslice vector scale

    int cAdd, cSub, cMul;       // conventional add/sub/mul
    int chAdd;                  // conventional horizontal adder
    int cvScal;                 // conventional vector scale
    int cDot;                   // conventional dot product

    int places = 2; //number of decimal places for fixed point function

    int count = 0;
    while (count < PASSES)
    {
        /* Test the conversion functionality with structures 'charlie' and 'delta' */
        std_type * charlie = (std_type *) calloc (sizeof(base_type)*8,sizeof(std_type));
        insert_random_values(charlie, 256, sizeof(base_type)*8);
        base_type * charlieBitslice = (base_type *) calloc (size,sizeof(base_type));
        start = __rdtsc();
        bitslice_convert(charlie,charlieBitslice,sizeof(base_type)*8);
        stop = __rdtsc();
        convert_timings[count] = stop - start;

        /*for (int i = 0; i < sizeof(base_type)*8; i++){
            std::cout << "charlie["<<i<<"]: " << (unsigned) charlie[i] << std::endl;
        }

        for (int i = 0; i < size; i++){
            std::cout << "charlieBitslice["<<i<<"]: " << std::bitset<sizeof(base_type)*8>(charlieBitslice[i]) << std::endl;
        }*/

        std_type * delta = (std_type *) calloc (sizeof(base_type)*8,sizeof(std_type));

        start = __rdtsc();
        bitslice_convert_back(charlieBitslice,delta,sizeof(std_type)*8);
        stop = __rdtsc();
        convert_back_timings[count] = stop - start;

        /*for (int i = 0; i < sizeof(base_type)*8; i++){
            std::cout << "delta["<<i<<"]: " << (unsigned) delta[i] << std::endl;
        }*/

        free(charlie);
        free(charlieBitslice);
        free(delta);


        /* Addition */
        start = __rdtsc();
        bitslice_add(a, b, sums, size);
        stop = __rdtsc();
        add_timings[count] = stop - start;

        /* Subtraction */
        start = __rdtsc();
        bitslice_sub(a, b, diffs, size);
        stop = __rdtsc();
        sub_timings[count] = stop - start;

        /* Multiplication */
        start = __rdtsc();
        bitslice_mul(a,b,products,size);
        stop = __rdtsc();
        mul_timings[count] = stop - start;

        /* Signed Multiplication */
        start = __rdtsc();
        bitslice_signed_mul(a,b,signed_result,size);
        stop = __rdtsc();
        signed_mul_timings[count] = stop - start;

        /* Fixed-Point Multiplication */
        start = __rdtsc();
        bitslice_fp_signed_mul(a,b,fixed_point_result,size,places);
        stop = __rdtsc();
        fixed_point_timings[count] = stop - start;

        /* Horizontal Add */
        start = __rdtsc();
        bitslice_horizontal_add(a, horizontal_add_result, size);
        stop = __rdtsc();
        horizontal_add_timings[count] = stop - start;

        /* Vector Scale */
        start = __rdtsc();
        bitslice_blas_scale(a,scale,scale_result,size);
        stop = __rdtsc();
        scale_timings[count] = stop - start;

        /* Dot product */
        start = __rdtsc();
        bitslice_blas_dot(a,b,dot_result,size);
        stop = __rdtsc();
        dot_timings[count] = stop - start;

        /* Normal addition */
        start = __rdtsc();
        for (int i = 0; i < sizeof(base_type)*8; i++){
            normal_add_results[i] = alpha[i] + beta[i];
        }
        stop = __rdtsc();
        normal_add_times[count] = stop - start;

        /* Normal subtraction */
        start = __rdtsc();
        for (int i = 0; i < sizeof(base_type)*8; i++){
            normal_sub_results[i] = alpha[i] - beta[i];
        }
        stop = __rdtsc();
        normal_sub_times[count] = stop - start;

        /* Normal multiplication */
        start = __rdtsc();
        for (int i = 0; i < sizeof(base_type)*8; i++){
            normal_mul_results[i] = alpha[i] * beta[i];
        }
        stop = __rdtsc();
        normal_mul_times[count] = stop - start;

        /* Normal horizontal addition */
        normal_hAdd_result = 0;
        start = __rdtsc();
        for (int i = 0; i < sizeof(base_type)*8; i++){
            normal_hAdd_result += alpha[i];
        }
        stop = __rdtsc();
        normal_hAdd_times[count] = stop - start;

        /* Normal vector scale */
        start = __rdtsc();
        for (int i = 0; i < sizeof(base_type)*8; i++){
            normal_vector_scale_result[i] = alpha[i] * 3;
        }
        stop = __rdtsc();
        normal_vector_scale_timings[count] = stop - start;

        /* Normal dot product */
        start = __rdtsc();
        for (int i = 0; i < sizeof(base_type)*8; i++){
            normal_dot_product_result += alpha[i] * beta[i];
        }
        stop = __rdtsc();
        normal_dot_product_timings[count] = stop - start;

        count++;
    }

    /* Print Bitslice Arithmetic Operations */
    std::cout << "\n\nUNSIGNED ADDITION (A + B = C)\n";
    print_bitslice_arrays_demo(a, b, sums, size);
    uAdd = (int) getMedian(add_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice addition: " << uAdd << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    print_clock_cycles(add_timings, PASSES);

    std::cout << "\n\nUNSIGNED SUBTRACTION (A - B = C)\n";
    print_bitslice_arrays_demo(a, b, diffs, size);
    uSub = (int) getMedian(sub_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice subtraction: " << uSub << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(sub_timings, PASSES);

    std::cout << "\n\nUNSIGNED MULTIPLICATION (A * B = C)\n";
    print_bitslice_arrays_demo(a, b, products, size);
    uMul = (int) getMedian(mul_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice unsigned multiplication: " << uMul << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(mul_timings, PASSES);

    std::cout << "\n\nSIGNED MULTIPLICATION RESULT of  (A * B):  Modified Baugh-Wooley Algorithm\n";
    print_bitslice_array(signed_result, size*2);
    sMul = (int) getMedian(signed_mul_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice signed multiplication: " << sMul << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(signed_mul_timings, PASSES);

    /* Fixed Point */
    std::cout << "\n\nFIXED-POINT MULTIPLICATION RESULT of (A * B)\n";
    print_bitslice_fp_array(fixed_point_result,size*2,places);
    fixp = (int) getMedian(fixed_point_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice fixed-point multiplication: " << fixp << std::endl;

    /* Bitslice Horizontal Adder */
    std::cout << "\n\nHORIZONTAL ADDITION\n";
    std::cout << "(the leftmost bitslice contains the sum of all bitslice values contained in A. The rest are garbage results)\n";
    print_bitslice_array(horizontal_add_result, size);
    hAdd = (int) getMedian(horizontal_add_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice horizontal addition: " << hAdd << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(horizontal_add_timings, PASSES);

    /* Bitslice BLAS: Vector Scale */
    std::cout << "\n\nBLAS: VECTOR SCALE (of A) - scaled by 3\n";
    print_bitslice_array(scale_result, size*2);
    scal = (int) getMedian(scale_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice vector scale: " << scal << std::endl;

    /* uncomment if you need individual clock cycles printing out */
   // print_clock_cycles(scale_timings, PASSES);

    /* Bitslice BLAS: Dot product */
    std::cout << "\n\nBLAS: DOT PRODUCT (of A and B)\n";
    std::cout << "(the leftmost bitslice contains the result. The rest are garbage results)\n";
    print_bitslice_array(dot_result, size*2);
    bDot = (int) getMedian(dot_timings, PASSES);
    std::cout << "\nMedian execution time for bitslice dot product: " << bDot << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(dot_timings, PASSES);

    /* Print Conventional Arithmetic Operations */
    std::cout << "\n\nCONVENTIONAL ADDITION\n";
    for (int i = 0; i < sizeof(base_type)*8; i++){
        std::cout << (int)normal_add_results[i] << ", ";
    }
    cAdd = (int) getMedian(normal_add_times, PASSES);
    std::cout << "\nMedian execution time for conventional addition: " << cAdd << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(normal_add_times, PASSES);

    std::cout << "\n\nCONVENTIONAL SUBTRACTION\n";
    for (int i = 0; i < sizeof(base_type)*8; i++){
        std::cout << (int)normal_sub_results[i] << ", ";
    }
    cSub = (int) getMedian(normal_sub_times, PASSES);
    std::cout << "\nMedian execution time for conventional subtraction: " << cSub << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(normal_sub_times, PASSES);

    std::cout << "\n\nCONVENTIONAL (SIGNED) MULTIPLICATION\n";
    for (int i = 0; i < sizeof(base_type)*8; i++){
        std::cout << (int)normal_mul_results[i] << ", ";
    }
    cMul = (int) getMedian(normal_mul_times, PASSES);
    std::cout << "\nMedian execution time for conventional multiplication: " << cMul << std::endl;

    /* uncomment if you need individual clock cycles printing out */
   // print_clock_cycles(normal_mul_times, PASSES);

    std::cout << "\n\nCONVENTIONAL HORIZONTAL ADDITION\n";
    std::cout << "Sum: " << (int)normal_hAdd_result << std::endl;
    chAdd = (int) getMedian(normal_hAdd_times, PASSES);
    std::cout << "Median execution time for conventional horizontal addition: " << chAdd << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(normal_hAdd_times, PASSES);

    std::cout << "\n\nCONVENTIONAL VECTOR SCALE\n";
    for (int i = 0; i < sizeof(base_type)*8; i++){
        std::cout << (int) normal_vector_scale_result[i] << ", ";
    }
    cvScal = (int) getMedian(normal_vector_scale_timings, PASSES);
    std::cout << "\nMedian execution time for conventional vector scale: " << cvScal << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(normal_vector_scale_timings, PASSES);

    std::cout << "\n\nCONVENTIONAL DOT PRODUCT\n";
    std::cout << "Dot product: " << normal_dot_product_result << std::endl;
    cDot = (int) getMedian(normal_dot_product_timings, PASSES);
    std::cout << "Median execution time for conventional dot product: " << cDot << std::endl;

    /* uncomment if you need individual clock cycles printing out */
    //print_clock_cycles(normal_vector_scale_timings, PASSES);

    /* conversion timings */
    int convert = (int) getMedian(convert_timings, PASSES);
    int convert_back = (int) getMedian(convert_back_timings, PASSES);

    /* Results Table */
    std::cout << "\n\nMedian Clock Cycles (after 100 passes);\n";
    std::cout << "        \tAdd\tSub\tMul\tsMul\tFixP\thAdd\tvScal\tDotP\n";
    std::cout << "bitslice\t" << uAdd << "\t" << uSub << "\t" << uMul << "\t" << sMul << "\t" << fixp << "\t" << hAdd << "\t" << scal << "\t" << bDot << std::endl;
    std::cout << "normal  \t" << cAdd << "\t" << cSub << "\t" << cMul << "\t" << cMul << "\t" << "n/a\t" << chAdd << "\t" << cvScal << "\t" << cDot << std::endl;

    std::cout << "\n\nMedian conversion timings;\n";

    /* uncomment if you need individual clock cycles printing out */
    //std::cout << "Convert To\n";
    //print_clock_cycles(convert_timings, PASSES);
    //std::cout << "\nConvert Back\n";
    //print_clock_cycles(convert_back_timings, PASSES);

    std::cout << "Conventional to Bitslice Conversion: " << convert << std::endl;
    std::cout << "Bitslice to Conventional Conversion: " << convert_back << std::endl;

    /* Clean up */
    free(sums);
    free(diffs);
    free(products);
    free(productsDouble);
    free(signed_result);
    free(fixed_point_result);
    free(normal_add_results);
    free(normal_sub_results);
    free(normal_mul_results);
    free(horizontal_add_result);
    free(dot_result);
    free(scale_result);
    free(normal_vector_scale_result);

    return 0;
}
