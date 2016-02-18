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



#define DIF_BUTTERFLY(xi,xj,n)                                \
  do {                                                  \
    v16 u= xi;					\
    v16 v= xj;					\
    xi =  v16_add(u, v);				\
    if (n)                                              \
      xj = v16_shift_l(v16_sub(u, v), 2*n);		\
    else                                                \
      xj = v16_sub(u, v);				\
  } while(0)

#define DIT_BUTTERFLY(xi,xj,n)                                \
  do {                                                  \
    v16 u= xi;                                        \
    v16 v= xj;                                        \
    if (n)                                              \
      v = v16_shift_l(v, 2*n);			\
    xi = v16_add(u, v);                               \
    xj = v16_sub(u, v);                               \
  } while(0)


#define INTERLEAVE(x,y) v16_interleave_inplace(x, y)

void dif_fft8(void *a) {
  v16* const A = a;
  register v16 X0, X1, X2, X3, X4, X5, X6, X7;

  X0 = A[0];
  X1 = A[1];
  X2 = A[2];
  X3 = A[3];
  X4 = A[4];
  X5 = A[5];
  X6 = A[6];
  X7 = A[7];

  // TODO : check bounds everywhere

  // outer layer
  DIF_BUTTERFLY(X0, X4, 0);
  DIF_BUTTERFLY(X1, X5, 1);
  DIF_BUTTERFLY(X2, X6, 2);
  DIF_BUTTERFLY(X3, X7, 3);
  
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);
  
  // middle layer
  DIF_BUTTERFLY(X0, X2, 0);
  DIF_BUTTERFLY(X4, X6, 0);
  DIF_BUTTERFLY(X1, X3, 2);
  DIF_BUTTERFLY(X5, X7, 2);
  
  // external layer
  DIF_BUTTERFLY(X0, X1, 0);
  DIF_BUTTERFLY(X2, X3, 0);
  DIF_BUTTERFLY(X4, X5, 0);
  DIF_BUTTERFLY(X6, X7, 0);
  
  /* We don't need to reduce X(0) */
  A[0] = REDUCE_FULL(X0); // FIXME later
  A[1] = REDUCE_FULL(X1);
  A[2] = REDUCE_FULL(X2);
  A[3] = REDUCE_FULL(X3);
  A[4] = REDUCE_FULL(X4);
  A[5] = REDUCE_FULL(X5);
  A[6] = REDUCE_FULL(X6);
  A[7] = REDUCE_FULL(X7);
}




/*
 * FFT64
 * Input in [-128,383]
 *   Input[0,1,2,4,5,6] can be in [-256,383]
 * Output in [-128,128]
 *
 */

void fft64(void *a) {

  v16* const A = a;

  register v16 X0, X1, X2, X3, X4, X5, X6, X7;

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

  //outer layer
  DIF_BUTTERFLY(X0, X4, 0);
  DIF_BUTTERFLY(X1, X5, 1);
  DIF_BUTTERFLY(X2, X6, 2);
  DIF_BUTTERFLY(X3, X7, 3);
  
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);
  
  // middle layer
  DIF_BUTTERFLY(X0, X2, 0);
  DIF_BUTTERFLY(X4, X6, 0);
  DIF_BUTTERFLY(X1, X3, 2);
  DIF_BUTTERFLY(X5, X7, 2);
  
  // external layer
  DIF_BUTTERFLY(X0, X1, 0);
  DIF_BUTTERFLY(X2, X3, 0);
  DIF_BUTTERFLY(X4, X5, 0);
  DIF_BUTTERFLY(X6, X7, 0);
 
  /* We don't need to reduce X(0) */
  X1 = REDUCE_FULL(X1);
  X2 = REDUCE_FULL(X2);
  X3 = REDUCE_FULL(X3);
  X4 = REDUCE_FULL(X4);
  X5 = REDUCE_FULL(X5);
  X6 = REDUCE_FULL(X6);
  X7 = REDUCE_FULL(X7);
    
  /*
   * Multiply by twiddle factors
   */

  X1 = v16_mul(X1, FFT64_Twiddle[0]);
  X2 = v16_mul(X2, FFT64_Twiddle[1]);
  X3 = v16_mul(X3, FFT64_Twiddle[2]);
  X4 = v16_mul(X4, FFT64_Twiddle[3]);
  X5 = v16_mul(X5, FFT64_Twiddle[4]);
  X6 = v16_mul(X6, FFT64_Twiddle[5]);
  X7 = v16_mul(X7, FFT64_Twiddle[6]);

  /*
   * Transpose the FFT state with a revbin order permutation
   * on the rows and the column.
   * This will make the full FFT_64 in order.
   */

  INTERLEAVE(X0, X1);
  INTERLEAVE(X2, X3);
  INTERLEAVE(X4, X5);
  INTERLEAVE(X6, X7);

  INTERLEAVE(X0, X2);
  INTERLEAVE(X1, X3);
  INTERLEAVE(X4, X6);
  INTERLEAVE(X5, X7);

  INTERLEAVE(X0, X4);
  INTERLEAVE(X1, X5);
  INTERLEAVE(X2, X6);
  INTERLEAVE(X3, X7);


  /*
   * Finish with 8 parallels DIT FFT_8
   *
   * FFT_8 using w=4 as 8th root of unity
   *  Unrolled decimation in time (DIT) radix-2 NTT.
   *  Intput data is in revbin_permuted order.
   */
  

  X0 = REDUCE(X0);
  X1 = REDUCE(X1);
  X2 = REDUCE(X2);
  X3 = REDUCE(X3);
  X4 = REDUCE(X4);
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);
  
  DIT_BUTTERFLY(X0, X1, 0);
  DIT_BUTTERFLY(X2, X3, 0);
  DIT_BUTTERFLY(X4, X5, 0);
  DIT_BUTTERFLY(X6, X7, 0);
  
  DIT_BUTTERFLY(X0, X2, 0);
  DIT_BUTTERFLY(X4, X6, 0);
  DIT_BUTTERFLY(X1, X3, 2);
  DIT_BUTTERFLY(X5, X7, 2);
  
  X7 = REDUCE(X7);
  
  DIT_BUTTERFLY(X0, X4, 0);
  DIT_BUTTERFLY(X1, X5, 1);
  DIT_BUTTERFLY(X2, X6, 2);
  DIT_BUTTERFLY(X3, X7, 3);
  
  A[0] = REDUCE_FULL(X0);
  A[1] = REDUCE_FULL(X1);
  A[2] = REDUCE_FULL(X2);
  A[3] = REDUCE_FULL(X3);
  A[4] = REDUCE_FULL(X4);
  A[5] = REDUCE_FULL(X5);
  A[6] = REDUCE_FULL(X6);
  A[7] = REDUCE_FULL(X7);
  

}

/*
 * FFT128
 * Input in [-128,383]
 * Outut in [-128,128]
 *
 */
static inline void fft128(void *a) {

  int i;

  // Temp space to help for interleaving in the end
  v16 B[8];

  v16 *A = (v16*) a;
  //  v16 *Twiddle = (v16*)FFT128_Twiddle;

  /* Size-2 butterflies */

  for (i = 0; i<8; i++) {
    B[i]   = v16_add(A[i], A[i+8]);
    A[i+8] = v16_sub(A[i], A[i+8]);
    A[i+8] = REDUCE_FULL(A[i+8]);
    A[i+8] = v16_mul(A[i+8], FFT128_Twiddle[i]);
    A[i+8] = REDUCE(A[i+8]);
  }

  B[3] = REDUCE(B[3]);
  B[7] = REDUCE(B[7]);
  fft64(B);
  fft64(A+8);

  /* Transpose (i.e. interleave) */
  v16 *A1=A+8, *B1=B;
  
  for (i=0; i<8; i++) {
    A[2*i]   = *(B1++);
    A[2*i+1] = *(A1++);
    v16_interleave_inplace(A[2*i],A[2*i+1]);
  }

}
