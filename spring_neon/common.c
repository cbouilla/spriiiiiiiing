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
const int NBITS = 4;
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

v16 rand_v16() {
  v16 x;
  for(int j=0; j < 8; j++) {
    do {
      x[j] = rand();
    } while(x[j]==0);
  }
  return REDUCE_FULL(x);
}

/*
 * Rouding functions.
 */

// On renvoie 4 bits par coefficients.
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

// méthode top-secrète pour initialiser A et les s_i. Ne pas divulguer au public ! (à remplacer plus tard par lfsr ???)
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


