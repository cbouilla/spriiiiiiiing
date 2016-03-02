#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.h"

/*
 * Code BCH
 */
uint64_t BCH128to64(uv64 in){
  register ui64 b1 = in[0];
  register ui64 res_a = b1;
  register ui64 res_b = vshl_n_u64(b1, 2);
  register ui64 res_c = vshl_n_u64(b1, 7);
  register ui64 res_d = vshl_n_u64(b1, 8);

  res_a = ui64_shiftl_xor(res_a, b1, 10);
  res_b = ui64_shiftl_xor(res_b, b1, 12);
  res_c = ui64_shiftl_xor(res_c, b1, 14);
  res_d = ui64_shiftl_xor(res_d, b1, 15);
  res_a = ui64_shiftl_xor(res_a, b1, 16);
  res_b = ui64_shiftl_xor(res_b, b1, 23);
  res_c = ui64_shiftl_xor(res_c, b1, 25);
  res_d = ui64_shiftl_xor(res_d, b1, 27);
  res_a = ui64_shiftl_xor(res_a, b1, 28);
  res_b = ui64_shiftl_xor(res_b, b1, 30);
  res_c = ui64_shiftl_xor(res_c, b1, 31);
  res_d = ui64_shiftl_xor(res_d, b1, 32);
  res_a = ui64_shiftl_xor(res_a, b1, 33);
  res_b = ui64_shiftl_xor(res_b, b1, 37);
  res_c = ui64_shiftl_xor(res_c, b1, 38);
  res_d = ui64_shiftl_xor(res_d, b1, 39);
  res_a = ui64_shiftl_xor(res_a, b1, 40);
  res_b = ui64_shiftl_xor(res_b, b1, 41);
  res_c = ui64_shiftl_xor(res_c, b1, 42);
  res_d = ui64_shiftl_xor(res_d, b1, 44);
  res_a = ui64_shiftl_xor(res_a, b1, 45);
  res_b = ui64_shiftl_xor(res_b, b1, 48);
  res_c = ui64_shiftl_xor(res_c, b1, 58);
  res_d = ui64_shiftl_xor(res_d, b1, 61);
  res_a = ui64_shiftl_xor(res_a, b1, 63);
  
  register ui64 b2 = in[1];
  res_b = ui64_shiftr_xor(res_b, b2, 62);
  res_c = ui64_shiftr_xor(res_c, b2, 57);
  res_d = ui64_shiftr_xor(res_d, b2, 56);
  res_a = ui64_shiftr_xor(res_a, b2, 54);
  res_b = ui64_shiftr_xor(res_b, b2, 52);
  res_c = ui64_shiftr_xor(res_c, b2, 50);
  res_d = ui64_shiftr_xor(res_d, b2, 49);
  res_a = ui64_shiftr_xor(res_a, b2, 48);
  res_b = ui64_shiftr_xor(res_b, b2, 41);
  res_c = ui64_shiftr_xor(res_c, b2, 39);
  res_d = ui64_shiftr_xor(res_d, b2, 37);
  res_a = ui64_shiftr_xor(res_a, b2, 36);
  res_b = ui64_shiftr_xor(res_b, b2, 34);
  res_c = ui64_shiftr_xor(res_c, b2, 33);
  res_d = ui64_shiftr_xor(res_d, b2, 32);
  res_a = ui64_shiftr_xor(res_a, b2, 31);
  res_b = ui64_shiftr_xor(res_b, b2, 27);
  res_c = ui64_shiftr_xor(res_c, b2, 26);
  res_d = ui64_shiftr_xor(res_d, b2, 25);
  res_a = ui64_shiftr_xor(res_a, b2, 24);
  res_b = ui64_shiftr_xor(res_b, b2, 23);
  res_c = ui64_shiftr_xor(res_c, b2, 22);
  res_d = ui64_shiftr_xor(res_d, b2, 20);
  res_a = ui64_shiftr_xor(res_a, b2, 19);
  res_b = ui64_shiftr_xor(res_b, b2, 16);
  res_c = ui64_shiftr_xor(res_c, b2, 6);
  res_d = ui64_shiftr_xor(res_d, b2, 3);
  res_a = ui64_shiftr_xor(res_a, b2, 1);

  ui64 res = res_a ^ res_b ^ res_c ^ res_d;

  return (uint64_t)res ^ ((uint64_t)(-(b2&1)));
}
