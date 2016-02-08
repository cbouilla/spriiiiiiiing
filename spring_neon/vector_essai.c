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

int main(){
  v16 a, b, c;
  int i;

  for(i = 0; i < 8; i++){
    a[i] = 1;
    b[i] = i;
  }

  c = v16_add(a,b);

  for(i = 0; i < 8; i++){
    printf("c[%d] = %d\n", i, c[i]);
  }

}

