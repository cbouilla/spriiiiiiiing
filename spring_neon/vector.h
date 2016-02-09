#ifndef __VECTOR_H__
#define __VECTOR_H__

#define CAT(x, y) x##y
#define XCAT(x,y) CAT(x,y)

#include <stdint.h>
#include <arm_neon.h>

typedef int16x8_t v16 __attribute__ ((aligned (16)));
typedef int8x16_t v8 __attribute__ ((aligned (16)));
typedef int8x8_t sv8 __attribute__ ((aligned (16)));
typedef int8x8x2_t dsv8 __attribute__ ((aligned (16)));
typedef int32x2_t sv32 __attribute__((aligned (16)));

#define v16_and vandq_s16
#define v16_or  vorrq_s16
#define v16_xor veorq_s16

#define v16_add vaddq_s16
#define v16_sub vsubq_s16
#define v16_mul vmulq_s16
#define v16_neg vnegq_s16

#define v16_shift_r  vshrq_n_s16
#define v16_shift_l  vshlq_n_s16

#define v16_to_v8 vreinterpretq_s8_s16
#define sv8_to_sv32 vreinterpret_s32_s8

#define v16_cmpt_gt(a,b) vcgtq_s16(a,b)
#define v16_cmpt_eq(a,b) vceqq_s16(a,b)

#define v8_lside vget_high_s8
#define v8_rside vget_low_s8

#define v8_transpose(x) vtrn_s8(v8_lside(x), v8_rside(x))

#define  v16_interleave_inplace(a__,b__) ({		\
      int16x8x2_t c__ = vzipq_s16 (a__, b__);		\
      a__ = c__.val[0];					\
      b__ = c__.val[1];					\
    })


#define  v16_merge_inplace(a__,b__) ({				\
      int16x8x2_t c__ = vzipq_s16 (V3216(a__), V3216(b__));	\
      a__ = V1632(c__.val[0]);					\
      b__ = V1632(c__.val[1]);					\
    })


#define v16_movemask(x) ({			\
      dsv8 q0 = v8_transpose(v16_to_v8(x));	\
      sv32 q1 = vreinterpret_s32_s8(q0.val[0]);	\
      q1[0]|q1[1];				\
    })
