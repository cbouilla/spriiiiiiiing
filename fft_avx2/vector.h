#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>
#include <immintrin.h>

typedef __v8hi v16;


#define v16_and      _mm_and_si128
#define v16_movemask _mm_movemask_epi8
#define v16_add      _mm_add_epi16
#define v16_sub      _mm_sub_epi16
#define v16_mul      _mm_mullo_epi16
#define v16_shift_l  _mm_slli_epi16
#define v16_shift_r  _mm_srai_epi16
#define v16_cmp      _mm_cmpgt_epi16
#define v16_cmp_eq   _mm_cmpeq_epi16

#define v16_interleavel   _mm_unpacklo_epi16
#define v16_interleaveh   _mm_unpackhi_epi16

#ifdef __AVX2INTRIN_H
  typedef __v16hi v32;
  #define v32_cst(x) {x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x}
  #define v32_cmp_gt   _mm256_cmpgt_epi16
  #define v32_cmp_eq   _mm256_cmpeq_epi16
  #define v32_movemask _mm256_movemask_epi8
  #define v32_shift_l  _mm256_slli_epi16
  #define v32_shift_r  _mm256_srai_epi16
#endif

// useless shit
#define CV(x) {{x, x, x, x, x, x, x, x}}
union cv {
  unsigned short u16[8];
  v16 v16;
};
#define CAT(x, y) x##y
#define XCAT(x,y) CAT(x,y)


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
