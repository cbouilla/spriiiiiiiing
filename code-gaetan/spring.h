#ifndef PRF_SWIFFT_SIMD_H_
#define PRF_SWIFFT_SIMD_H_


#include <inttypes.h>

// K can be chosen from {64,128}.
// N must be 128
#define N 128
#define K 128

// Align to cache line size
#define ALIGN __attribute__ ((aligned (64)))


extern const uint8_t _randGeneratedNumbersForRq[] ALIGN;
extern const uint8_t _randGeneratedNumbersForR2[] ALIGN;

// number of iterations to run the method. For benchmarking
//#define PRINT_OUTPUT 1
//#define N_ITERATIONS 100000000
#define N_ITERATIONS 40000000
//#define N_ITERATIONS 8

/***** Preformance tweak *****/

/** Use CLMUL for polynomial multiplication [Autodetected] **/
// #define USE_CLMUL 1

/** Use CLMUL instead of subset-sum in R2 [Autodetected] **/
// #define USE_NAIVE_R2 1

/** Use precomputed tables for subset-sum [Should be 1] **/
#define SUBSET_TABLES 1

/** Use prefeteching to get key data in cache faster [Should be 1] **/
#define PREFETCH 1

/** Use precomputed tables for first layer of ExpToRadix **/
#define EXP_TABLES 1

// Autodetect CLMUL.
#ifndef USE_CLMUL
#ifdef __PCLMUL__
#define USE_CLMUL 1
#else
#define USE_CLMUL 0
#endif
#endif

#ifndef USE_NAIVE_R2
#define USE_NAIVE_R2 USE_CLMUL
#endif

#if USE_NAIVE_R2 && !USE_CLMUL
#warning USE_NAIVE_R2 is VERY slow without PCLMUL support
#endif

#if USE_NAIVE_R2 && !SUBSET_TABLES
#warning USE_NAIVE_R2 needs SUBSET_TABLES
#endif

#endif /* PRF_SWIFFT_SIMD_H_ */
