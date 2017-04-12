#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.c"

v16 A[16] __attribute__ ((aligned(16)));
v16 S_Eval[64][2][16] __attribute__ ((aligned(16)));

// a table of the powers of omega. 
const v16 omegaPowers[16] = {
    {1, -94, 98, 40, 95, 65, 58, -55},
    {30, 7, 113, -85, 23, -106, -59, -108},
    {-128, -47, 49, 20, -81, -96, 29, 101},
    {15, -125, -72, 86, -117, -53, 99, -54},
    {-64, 105, -104, 10, 88, -48, -114, -78},
    {-121, 66, -36, 43, 70, 102, -79, -27},
    {-32, -76, -52, 5, 44, -24, -57, -39},
    {68, 33, -18, -107, 35, 51, 89, 115},
    {-16, -38, -26, -126, 22, -12, 100, 109},
    {34, -112, -9, 75, -111, -103, -84, -71},
    {-8, -19, -13, -63, 11, -6, 50, -74},
    {17, -56, 124, -91, 73, 77, -42, 93},
    {-4, 119, 122, 97, -123, -3, 25, -37},
    {-120, -28, 62, 83, -92, -90, -21, -82},
    {-2, -69, 61, -80, 67, 127, -116, 110},
    {-60, -14, 31, -87, -46, -45, 118, -41}
};

const v16 REJECTION_VECT = CV(-1);
const v16 VECT_ZERO = CV(0);
const int NBITS = 4;
const v16 m = {0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0};
const v16 m0 = {0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00};
const v16 m1 = {0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0};
const sv16 sm0 = {0xff, 0xff, 0x00, 0x00};
const sv16 sm1 = {0x00, 0x00, 0xff, 0xff};
const v16 tm0 = {0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00};
const v16 tm1 = {0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0};

// Get Coefficients from FFT.
void ConvertEvalToCoefficients(const v16 *Eval, v16 *Coef) {
  for(int i = 0; i < 16; i++){
    Coef[i] = Eval[i];
  }

  fft128(Coef);

  for(int i = 0; i < 16; i++){
    Coef[i] = REDUCE_FULL(Coef[i] * omegaPowers[i]);
  }
}

uint32_t UpdateGray(uint32_t x, uint32_t Gray){

int flip = __builtin_ctz(Gray);
 int32_t mask;

  mask = (1 << flip);
  x ^= mask;

  return x;
}


/*
 * Update the Gray counter and multiply polynomials.
 */
uint32_t UpdateCounterMode(uint32_t x, v16 *Prod, uint32_t Gray){
  int flip;
  uint32_t mask, inv;

  //Index of the flipped bit
  flip = __builtin_ctz(Gray);

  mask = (1 << flip);
  x ^= mask;

  // New value of the flipped bit
  inv = (x >> flip) & 1;

  // Which s_i to use.

  const v16 *s = S_Eval[flip][inv];

  for(int i = 0; i < 16; i++){
    Prod[i] = REDUCE_FULL(v16_mul(Prod[i], s[i]));
  }

  return x;
}

//rejection sampling
int reject(v16 a){
  v16 mask = v16_cmp_eq(a, REJECTION_VECT);
  return v16_movemask(mask);
}

/*
 * Return MSB of each coefficient of the mask (output is permuted).
 */
char PermutedMovmask16(v16 a){ 
  v8 b = vnegq_s8(vreinterpretq_s8_s16(a));
  dsv8 c = vzip_s8(vget_low_s8(b), vget_high_s8(b));
  sv8 d = vsli_n_s8(c.val[0], c.val[1], 4);
  sv32 e = vreinterpret_s32_s8(d);
  int32_t r0 = e[0] ^ (e[1] << 2);
  int16_t r1 = (r0 >>15) & (0x00ff);
  char r = r1 ^ (r0 >> 24);
  return r;
 }

/*
 * Extract most signifiquant bit of each coefficient.
 * Output is permuted.
 */
char PermutedMSB(v16 a){
  v16 mask = v16_cmp_gt(VECT_ZERO, a);
 char  r = PermutedMovmask16(mask);
  return r;
}


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


/*
 * Rounding functions.
 */

// Return 4 bits per coefficients (no permutation).
uint32_t rounding4(v16 a) {
  v16 b = v16_and(a, m0);
  v16 c = v16_shift_r(v16_and(a, m1), 4);
  dsv16 d0 = v16_transpose(b);
  dsv16 d1 = v16_transpose(c);
  sv16 d = sv16_xor(d0.val[0], d1.val[1]);
  sv16 e = sv16_shift_l(sv16_and(d, sm0), 8);
  sv16 f = sv16_and(d, sm1);
  uint32_t res = ((e[1]^f[3])<<16);
  res ^=(e[0]^f[2]); 

  return res;
}

// Return 4 bits per coefficients (permuted output).
uint32_t prounding4(v16 a){
  v8 b = vreinterpretq_s8_s16(a);
  dsv8 c = vzip_s8(vget_low_s8(b), vget_high_s8(b));
  sv8 d = vsri_n_s8(c.val[0], c.val[1], 4);
  sv32 e = vreinterpret_s32_s8(d);
  uint32_t r = (e[0]<<16)^(e[1]);
  return r;
}

uint16_t rounding2(v16 a) {
  v16 b = v16_and(a, tm0);
  v16 c = v16_shift_r(v16_and(a, tm1), 2);
  dsv16 d0 = v16_transpose(b);
  dsv16 d1 = v16_transpose(c);
  sv16 d = sv16_xor(d0.val[0], d1.val[1]);
  uint16_t res = (d[1]<<8)^(d[3]<<4);
  res ^= (d[0])^(d[2]>>4);
  return res;
}

v16 rand_v16() {
  v16 x;
  for(int j=0; j < 8; j++) {
    do {
      x[j] = rand();
    } while(x[j]==0);
  }
  return REDUCE_FULL(x);
}


//For benchmark. Do not use this for cryptographic purpose.
void init_secrets() {
  srand(42);

  for(int i=0; i < 64; i++) {
    for(int k=0; k<2; k++){
      for(int j=0; j < 16; j++) {
        S_Eval[i][k][j] = rand_v16();
        // OK, Sinv isn't really S inverse. And So ?
      }
    }
  }
  for(int i=0; i<16; i++) {
    A[i]=rand_v16();
  }
}


