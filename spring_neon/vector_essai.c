#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <arm_neon.h>

typedef int16x8_t v16 __attribute__ ((aligned (16)));

#define v16_and vandq_s16
#define v16_or  vorrq_s16
#define v16_xor veorq_s16

#define v16_add vaddq_s16
#define v16_sub vsubq_s16
#define v16_mul vmulq_s16

#define v16_shift_r  vshrq_n_s16
#define v16_shift_l  vshlq_n_s16

#define v16_movemask(x) ({						\
      int8x16_t q0 = vreinterpretq_s8_s16(a);				\
      int8x8x2_t q1 = vtrn_s8(vget_high_s8(q0), vget_low_s8(q0));	\
      int32x2_t d = vreinterpret_s32_s8(q1.val[1]);			\
      d[0]|d[1];							\
    })


int main(){
  v16 a;
  int i, mask;

  for(i = 0; i < 8; i++){
    a[i] = 0;
  }

  a[0] = 1;
  a[7] = 1;

  //c = v16_add(a,b);

  mask = v16_movemask(a);

  printf("mask = %x", mask);
  /* for(i = 0; i < 8; i++){ */
  /*   printf("c[%d] = %d\n", i, c[i]); */
  /* } */

}

