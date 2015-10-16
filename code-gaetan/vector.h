#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "compat.h"
#define CAT(x, y) x##y
#define XCAT(x,y) CAT(x,y)

// and shift_r shift_l sub cmp add mul v16_interleavel v16_interleaveh

#include <emmintrin.h>

typedef __v8hi v16;

#define CV(x) {{x, x, x, x, x, x, x, x}}

#define v16_and _mm_and_si128

#define v16_add      _mm_add_epi16
#define v16_sub      _mm_sub_epi16
#define v16_mul      _mm_mullo_epi16
#define v16_shift_l  _mm_slli_epi16
#define v16_shift_r  _mm_srai_epi16
#define v16_cmp      _mm_cmpgt_epi16

#define v16_interleavel   _mm_unpacklo_epi16
#define v16_interleaveh   _mm_unpackhi_epi16

/* Unions to convert vector types to scalar types
 */

union cv {
  unsigned short u16[8];
  v16 v16;
};

#endif // #ifndef __VECTOR_H__
