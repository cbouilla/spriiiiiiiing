#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>
#include <immintrin.h>
#include <wmmintrin.h>

typedef __v8hi v16;
typedef __v16hi v32;

#define v32_cst(x) {x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x}
/*#define v32_cmp_gt   _mm256_cmpgt_epi16
  #define v32_cmp_eq   _mm256_cmpeq_epi16
  #define v32_movemask _mm256_movemask_epi8
  #define v32_shift_l  _mm256_slli_epi16
  #define v32_shift_r  _mm256_srai_epi16*/

/*
 * Reduce modulo 257; result is in [-127; 383]
 * REDUCE(x) := (x&255) - (x>>8)
 */
#define REDUCE(x)                               \
  ((x & cst255) - _mm256_srai_epi16(x, 8))

/*
 * Reduce from [-127; 383] to [-128; 128]
 * EXTRA_REDUCE_S(x) := x<=128 ? x : x-257
 */
#define EXTRA_REDUCE(x)                       \
  (x - (cst257 & _mm256_cmpgt_epi16(x, cst128)))

/*
 * Reduce modulo 257; result is in [-128; 128]
 */
#define REDUCE_FULL(x)                        \
  EXTRA_REDUCE(REDUCE(x))


// useless shit
/*#define CV(x) {{x, x, x, x, x, x, x, x}}
union cv {
  unsigned short u16[8];
  v16 v16;
};
#define CAT(x, y) x##y
#define XCAT(x,y) CAT(x,y)
*/

static const v32 cst128 = v32_cst(128);
static const v32 cst255 = v32_cst(255);
static const v32 cst257 = v32_cst(257);


extern void fft128(void *a);
extern int reject(const v32 a);
extern uint16_t msb(const v32 a);
extern uint64_t rounding(const v32 a);
extern uint64_t BCH128to64_clmul (const __v2di in);

#define rdtsc()                                                         \
  ({                                                                    \
    uint32_t lo, hi;                                                    \
    __asm__ __volatile__ (      /* serialize */                         \
                          "xorl %%eax,%%eax \n        cpuid"            \
                          ::: "%rax", "%rbx", "%rcx", "%rdx");          \
    /* We cannot use "=A", since this would use %rax on x86_64 */       \
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));              \
    (uint64_t) hi << 32 | lo;                                            \
  })                                   

#endif // #ifndef __VECTOR_H__
