#include <stdlib.h>
#include <stdio.h>

#include "vector.h"

/* Twiddle tables */

  static const v16 FFT64_Twiddle[] = {
    {1,    2,    4,    8,   16,   32,   64,  128},
    {1,   60,    2,  120,    4,  -17,    8,  -34},
    {1,  120,    8,  -68,   64,  -30,   -2,   17},
    {1,   46,   60,  -67,    2,   92,  120,  123},
    {1,   92,  -17,  -22,   32,  117,  -30,   67},
    {1,  -67,  120,  -73,    8,  -22,  -68,  -70},
    {1,  123,  -34,  -70,  128,   67,   17,   35},
  };


  static const v16 FFT128_Twiddle[] =  {
    {  1, -118,   46,  -31,   60,  116,  -67,  -61},
    {  2,   21,   92,  -62,  120,  -25,  123, -122},
    {  4,   42,  -73, -124,  -17,  -50,  -11,   13},
    {  8,   84,  111,    9,  -34, -100,  -22,   26},
    { 16,  -89,  -35,   18,  -68,   57,  -44,   52},
    { 32,   79,  -70,   36,  121,  114,  -88,  104},
    { 64,  -99,  117,   72,  -15,  -29,   81,  -49},
    {128,   59,  -23, -113,  -30,  -58,  -95,  -98},
  };


  static const v16 FFT256_Twiddle[] =  {
    {   1,   41, -118,   45,   46,   87,  -31,   14}, 
    {  60, -110,  116, -127,  -67,   80,  -61,   69}, 
    {   2,   82,   21,   90,   92,  -83,  -62,   28}, 
    { 120,   37,  -25,    3,  123,  -97, -122, -119}, 
    {   4,  -93,   42,  -77,  -73,   91, -124,   56}, 
    { -17,   74,  -50,    6,  -11,   63,   13,   19}, 
    {   8,   71,   84,  103,  111,  -75,    9,  112}, 
    { -34, -109, -100,   12,  -22,  126,   26,   38}, 
    {  16, -115,  -89,  -51,  -35,  107,   18,  -33}, 
    { -68,   39,   57,   24,  -44,   -5,   52,   76}, 
    {  32,   27,   79, -102,  -70,  -43,   36,  -66}, 
    { 121,   78,  114,   48,  -88,  -10,  104, -105}, 
    {  64,   54,  -99,   53,  117,  -86,   72,  125}, 
    { -15, -101,  -29,   96,   81,  -20,  -49,   47}, 
    { 128,  108,   59,  106,  -23,   85, -113,   -7}, 
    { -30,   55,  -58,  -65,  -95,  -40,  -98,   94}
  };


/* static inline void fft128(void *a); */
 void fft64(void *a);


/*
 * Reduce modulo 257; result is in [-127; 383]
 * REDUCE(x) := (x&255) - (x>>8)
 */
#define REDUCE(x)                               \
  v16_sub(v16_and(x, V255), v16_shift_r (x, 8))

/*
 * Reduce from [-127; 383] to [-128; 128]
 * EXTRA_REDUCE_S(x) := x<=128 ? x : x-257
 */
#define EXTRA_REDUCE(x)                       \
  v16_sub(x, v16_and(V257, v16_cmp_gt(x, V128)))

/*
 * Reduce modulo 257; result is in [-128; 128]
 */
#define REDUCE_FULL(x)                        \
  EXTRA_REDUCE(REDUCE(x))


#define DO_REDUCE(i)                            \
  X(i) = REDUCE(X(i))

#define DO_REDUCE_FULL_S(i)                     \
  do {                                          \
    X(i) = REDUCE(X(i));                        \
    X(i) = EXTRA_REDUCE(X(i));                \
  } while(0)


void dif_fft8(void *a) {
  v16* const A = a;
  register v16 X0, X1, X2, X3, X4, X5, X6, X7;

#define X(i) X##i

  X0 = A[0];
  X1 = A[1];
  X2 = A[2];
  X3 = A[3];
  X4 = A[4];
  X5 = A[5];
  X6 = A[6];
  X7 = A[7];


  /*
   * Begin with 8 parallels DIF FFT_8
   *
   * FFT_8 using w=4 as 8th root of unity
   *  Unrolled decimation in frequency (DIF) radix-2 NTT.
   *  Output data is in revbin_permuted order.
   */

  #define wn0 0
  #define wn1 2
  #define wn2 4
  #define wn3 6

#define BUTTERFLY(i,j,n)                                \
  do {                                                  \
    v16 u= X(i);                                        \
    v16 v= X(j);                                        \
    X(i) =  v16_add(u, v);                              \
    if (n)                                              \
      X(j) = v16_shift_l(v16_sub(u, v), XCAT(wn,n));  \
    else                                                \
      X(j) = v16_sub(u, v);                             \
  } while(0)

  // outer layer
  BUTTERFLY(0, 4, 0);
  BUTTERFLY(1, 5, 1);
  BUTTERFLY(2, 6, 2);
  BUTTERFLY(3, 7, 3);
  
  DO_REDUCE(5);
  DO_REDUCE(6);
  DO_REDUCE(7);
  
  // middle layer
  BUTTERFLY(0, 2, 0);
  BUTTERFLY(4, 6, 0);
  BUTTERFLY(1, 3, 2);
  BUTTERFLY(5, 7, 2);
  
  // external layer
  BUTTERFLY(0, 1, 0);
  BUTTERFLY(2, 3, 0);
  BUTTERFLY(4, 5, 0);
  BUTTERFLY(6, 7, 0);
  
  /* We don't need to reduce X(0) */
  DO_REDUCE_FULL_S(1);
  DO_REDUCE_FULL_S(2);
  DO_REDUCE_FULL_S(3);
  DO_REDUCE_FULL_S(4);
  DO_REDUCE_FULL_S(5);
  DO_REDUCE_FULL_S(6);
  DO_REDUCE_FULL_S(7);
  
  A[0] = X0;
  A[1] = X1;
  A[2] = X2;
  A[3] = X3;
  A[4] = X4;
  A[5] = X5;
  A[6] = X6;
  A[7] = X7;

  #undef BUTTERFLY
}
