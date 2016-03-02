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

// Coef <--- coefficient du polynome représenté sous forme évalué dans Eval.
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
 * Met à jour le "Gray counter" et effectue la multiplication des évaluations.
 * x : la chaine de bits représentant la position du compteur.
 * Gray : le numéro de la valeur du compteur + 1 (déjà mis à jour)
 */
uint32_t UpdateCounterMode(uint32_t x, v16 *Prod, uint32_t Gray){
  int flip;
  uint32_t mask, inv;

  //index du bit qui doit changer
  flip = __builtin_ctz(Gray);

  mask = (1 << flip);
  x ^= mask;

  // nouvelle valeur du bit qui a changé
  inv = (x >> flip) & 1;

  // déterminer quel si à utiliser.

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

// Sans utiliser NEON


/*
 * renvoie le bit de poids fort de chaque coefficient d'un mask
 * la sortie est permutée.
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

/*
 * Code BCH
 */
uint64_t BCH128to64(uv64 in){
  register ui64 b1 = in[0];
  register ui64 res = b1;
  res = ui64_shiftl_xor(res, b1, 2);
  res = ui64_shiftl_xor(res, b1, 7);
  res = ui64_shiftl_xor(res, b1, 8);
  res = ui64_shiftl_xor(res, b1, 10);
  res = ui64_shiftl_xor(res, b1, 12);
  res = ui64_shiftl_xor(res, b1, 14);
  res = ui64_shiftl_xor(res, b1, 15);
  res = ui64_shiftl_xor(res, b1, 16);
  res = ui64_shiftl_xor(res, b1, 23);
  res = ui64_shiftl_xor(res, b1, 25);
  res = ui64_shiftl_xor(res, b1, 27);
  res = ui64_shiftl_xor(res, b1, 28);
  res = ui64_shiftl_xor(res, b1, 30);
  res = ui64_shiftl_xor(res, b1, 31);
  res = ui64_shiftl_xor(res, b1, 32);
  res = ui64_shiftl_xor(res, b1, 33);
  res = ui64_shiftl_xor(res, b1, 37);
  res = ui64_shiftl_xor(res, b1, 38);
  res = ui64_shiftl_xor(res, b1, 39);
  res = ui64_shiftl_xor(res, b1, 40);
  res = ui64_shiftl_xor(res, b1, 41);
  res = ui64_shiftl_xor(res, b1, 42);
  res = ui64_shiftl_xor(res, b1, 44);
  res = ui64_shiftl_xor(res, b1, 45);
  res = ui64_shiftl_xor(res, b1, 48);
  res = ui64_shiftl_xor(res, b1, 58);
  res = ui64_shiftl_xor(res, b1, 61);
  res = ui64_shiftl_xor(res, b1, 63);

  register ui64 b2 = in[1];
  res = ui64_shiftr_xor(res, b2, 62);
  res = ui64_shiftr_xor(res, b2, 57);
  res = ui64_shiftr_xor(res, b2, 56);
  res = ui64_shiftr_xor(res, b2, 54);
  res = ui64_shiftr_xor(res, b2, 52);
  res = ui64_shiftr_xor(res, b2, 50);
  res = ui64_shiftr_xor(res, b2, 49);
  res = ui64_shiftr_xor(res, b2, 48);
  res = ui64_shiftr_xor(res, b2, 41);
  res = ui64_shiftr_xor(res, b2, 39);
  res = ui64_shiftr_xor(res, b2, 37);
  res = ui64_shiftr_xor(res, b2, 36);
  res = ui64_shiftr_xor(res, b2, 34);
  res = ui64_shiftr_xor(res, b2, 33);
  res = ui64_shiftr_xor(res, b2, 32);
  res = ui64_shiftr_xor(res, b2, 31);
  res = ui64_shiftr_xor(res, b2, 27);
  res = ui64_shiftr_xor(res, b2, 26);
  res = ui64_shiftr_xor(res, b2, 25);
  res = ui64_shiftr_xor(res, b2, 24);
  res = ui64_shiftr_xor(res, b2, 23);
  res = ui64_shiftr_xor(res, b2, 22);
  res = ui64_shiftr_xor(res, b2, 20);
  res = ui64_shiftr_xor(res, b2, 19);
  res = ui64_shiftr_xor(res, b2, 16);
  res = ui64_shiftr_xor(res, b2, 6);
  res = ui64_shiftr_xor(res, b2, 3);
  res = ui64_shiftr_xor(res, b2, 1);

  return (uint64_t)res ^ ((uint64_t)(-(b2&1)));
}

/*
 * Rouding functions.
 */

// On renvoie 4 bits par coefficients (pas de permutations).
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

// on renvoie 4 bits par coefficients (sortie permutée).
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

/* uint32_t RandomSequence(uint32_t *X){ */
/*   uint32_t r; */

/*   r = X[3]; */
/*   X[3] = X[2]; */
/*   X[2] = X[1]; */
/*   X[1] = X[0]; */
/*   X[0] = X[0] * 0xdeadbeef; */
/*   X[0] += r; */

/*   return r;   */
/* } */

/* v16 Randomv16AndUpdateSequence(uint32_t *seed){ */
/*   v16 x; */
/*   int16_t x0, x1, a0, a1; */
/*   for(int i = 0; i < 7; i+=2){ */
/*     do{ */
/*       int32_t xx = RandomSequence(seed); */
/*       x0 = xx; x1 = xx>>16; */
/*       a0 = x0/256; a1 = x1/256; */
/*     } while(256*a0 == x0 || 256*a1 == x1); */
/*     x[i] = x0; */
/*     x[i+1] = x1; */
/*   } */
/*   return REDUCE_FULL(x); */
/* } */


// méthode top-secrète pour initialiser A et les s_i. Ne pas divulguer au public !
void init_secrets() {
  srand(42);

  for(int i=0; i < 64; i++) {
    for(int k=0; k<2; k++){
      for(int j=0; j < 16; j++) {
        S_Eval[i][k][j] = rand_v16();
        // OK, Sinv n'est pas vraiment l'inverse de S. Et alors ?
      }
    }
  }
  for(int i=0; i<16; i++) {
    A[i]=rand_v16();
  }
}


