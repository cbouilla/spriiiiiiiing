#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.h"
#include "vector32.c"

#define EXTRACTED_BITS 4
#define N_BYTES 320000000


v32 A[8];
v32 S_Eval[64][2][8];


// a table of the powers of omega. 
const v32 omegaPowers[8] = {
    {1, -94, 98, 40, 95, 65, 58, -55, 30, 7, 113, -85, 23, -106, -59, -108},
    {-128, -47, 49, 20, -81, -96, 29, 101, 15, -125, -72, 86, -117, -53, 99, -54},
    {-64, 105, -104, 10, 88, -48, -114, -78, -121, 66, -36, 43, 70, 102, -79, -27},
    {-32, -76, -52, 5, 44, -24, -57, -39, 68, 33, -18, -107, 35, 51, 89, 115},
    {-16, -38, -26, -126, 22, -12, 100, 109, 34, -112, -9, 75, -111, -103, -84, -71},
    {-8, -19, -13, -63, 11, -6, 50, -74, 17, -56, 124, -91, 73, 77, -42, 93},
    {-4, 119, 122, 97, -123, -3, 25, -37, -120, -28, 62, 83, -92, -90, -21, -82},
    {-2, -69, 61, -80, 67, 127, -116, 110, -60, -14, 31, -87, -46, -45, 118, -41}
};


// Coef <--- coefficients du polynôme représenté sous forme évalué dans Eval
void ConvertEvalToCoefficients(const v32 *Eval, v32 *Coef) {
  for(int i=0; i<8; i++) {
    Coef[i] = Eval[i];
  }

  fft128(Coef);

  for(int i=0; i<8; i++){
    Coef[i] = REDUCE_FULL(Coef[i] * omegaPowers[i]);
  }
}

/**
 * @param x : la chaine de bits représentant la position actuelle du "gray counter"
 * @param Gray : le numéro de la valeur générée
 */
uint64_t UpdateCounterMode(uint64_t x, v32 *Prod, const uint64_t Gray){
  int flip;
  uint64_t mask, inv;


  // index du bit qui doit changer
  flip = 1 + __builtin_ctz(Gray);

  mask = (1 << (flip));
  x ^= mask;

  // nouvelle valeur du bit qui a changé
  inv = (x >> flip) & 1;

  // jeu de s_i à utiliser
  //const v16 *s = inv ? Sinv_Eval[flip] : S_Eval[flip];
  const v32 *s = S_Eval[flip][inv];

  for(int i=0; i<8; i++) {
    Prod[i] = REDUCE_FULL(Prod[i] * s[i]);
  }
  return x;
}



int16_t nice_rand() {
  int x = rand() % 257;  
  return (x > 128) ? x - 257 : x;
}

// méthode top-secrète pour initialiser A et les s_i. Ne pas divulguer au public !
void init_secrets() {
  srand(42);

  for(int i=0; i < 64; i++) {
    for(int k=0; k<2; k++){
      int16_t * s_ = (int16_t *) & S_Eval[i][k][0];

      for(int j=0; j < 128; j++) {
        s_[j] = nice_rand();
      }
    }
  }

  int16_t * a_ = (int16_t *) &A[0];
  for(int i=0; i<128; i++) {
    a_[i] = nice_rand();
  }
}

/**************************** CODE DEDICATED TO THE PRFs ****************************/
vLog A_log[4];
vLog S_log[64][4];
vLog SS_tables[8][256][4];

void init_subset_sum_tables() {

  for(int i=0; i<8; i++) { // for each 8-bit pack in x
    for(int x=0; x<256; x++) { // for each possible 8-bit value

      // include A_log
      SS_tables[i][x][0] = i ? ZERO_VECT : A_log[0];
      SS_tables[i][x][1] = i ? ZERO_VECT : A_log[1];
      SS_tables[i][x][2] = i ? ZERO_VECT : A_log[2];
      SS_tables[i][x][3] = i ? ZERO_VECT : A_log[3];

      for(int k=0; k<8; k++) {
        if (x & (1 << k)) {
          SS_tables[i][x][0] += S_log[8*i + k][0];
          SS_tables[i][x][1] += S_log[8*i + k][1];
          SS_tables[i][x][2] += S_log[8*i + k][2];
          SS_tables[i][x][3] += S_log[8*i + k][3];
        }
      }
    }
  }
}


// a table of the powers of GENERATOR.  i --> GENERATOR ** I   ( mod FIELD_SIZE )
static const int16_t generatorPowers[257] = {
  1, 3, 9, 27, 81, -14, -42, -126, -121, -106, -61, 74, -35, -105,
  -58, 83, -8, -24, -72, 41, 123, 112, 79, -20, -60, 77, -26, -78,
  23, 69, -50, 107, 64, -65, 62, -71, 44, -125, -118, -97, -34,
  -102, -49, 110, 73, -38, -114, -85, 2, 6, 18, 54, -95, -28, -84,
  5, 15, 45, -122, -109, -70, 47, -116, -91, -16, -48, 113, 82,
  -11, -33, -99, -40, -120, -103, -52, 101, 46, -119, -100, -43,
  128, 127, 124, 115, 88, 7, 21, 63, -68, 53, -98, -37, -111, -76,
  29, 87, 4, 12, 36, 108, 67, -56, 89, 10, 30, 90, 13, 39, 117, 94,
  25, 75, -32, -96, -31, -93, -22, -66, 59, -80, 17, 51, -104, -55,
  92, 19, 57, -86, -1, -3, -9, -27, -81, 14, 42, 126, 121, 106, 61,
  -74, 35, 105, 58, -83, 8, 24, 72, -41, -123, -112, -79, 20, 60,
  -77, 26, 78, -23, -69, 50, -107, -64, 65, -62, 71, -44, 125, 118,
  97, 34, 102, 49, -110, -73, 38, 114, 85, -2, -6, -18, -54, 95,
  28, 84, -5, -15, -45, 122, 109, 70, -47, 116, 91, 16, 48, -113,
  -82, 11, 33, 99, 40, 120, 103, 52, -101, -46, 119, 100, 43, -128,
  -127, -124, -115, -88, -7, -21, -63, 68, -53, 98, 37, 111, 76,
  -29, -87, -4, -12, -36, -108, -67, 56, -89, -10, -30, -90, -13,
  -39, -117, -94, -25, -75, 32, 96, 31, 93, 22, 66, -59, 80, -17,
  -51, 104, 55, -92, -19, -57, 86, 1};


// Coef <--- coefficients du polynôme représenté sous forme de log d'évaluations dans Log
void ConvertSubsetSumToCoefficients(const vLog *Log, v32 *Coef) {
  
  // remplace les logs par les "vraies" valeurs. Eventuellement parallelisable.
  uint8_t *l_ = (uint8_t *) Log;
  int16_t *c_ = (int16_t *) Coef;
  for(int i = 0; i < 128; i++) {
    c_[i] = generatorPowers[ l_[i] ];
  }

  fft128(Coef);

  for(int i=0; i<8; i++){
    Coef[i] = REDUCE_FULL(Coef[i] * omegaPowers[i]);
  }
}



void inline ComputeSubsetSum(uint64_t x, vLog *Log) {
  // initialize Log with Log_a
    Log[0] = A_log[0];
    Log[1] = A_log[1];
    Log[2] = A_log[2];
    Log[3] = A_log[3];

#define likely(x)    __builtin_expect (!!(x), 1)

  // scan the bits of x, least signifiant to most significant
  while(likely(x)) {
    // clear the rightmost bit of x
    int j = __builtin_ffsll(x) - 1;
    x ^= (1ULL << j);

    Log[0] += S_log[j][0];
    Log[1] += S_log[j][1];
    Log[2] += S_log[j][2];
    Log[3] += S_log[j][3];
 
  }
}

void  ComputeSubsetSum_tabulated(uint64_t x, vLog *Log) {
  // initialize Log with Log_a
    
  int k = (x & 0xff);
  Log[0] = SS_tables[0][k][0];
  Log[1] = SS_tables[0][k][1];
  Log[2] = SS_tables[0][k][2];
  Log[3] = SS_tables[0][k][3];

  for(int i=1; i<8; i++) {
    k = (x >> (8ULL*i)) & 0xff;
    
    Log[0] += SS_tables[i][k][0];
    Log[1] += SS_tables[i][k][1];
    Log[2] += SS_tables[i][k][2];
    Log[3] += SS_tables[i][k][3];
  }
}


// méthode top-secrète pour initialiser A et les s_i. Ne pas divulguer au public !
void init_secrets_log() {
  srand(42);

  for(int i=0; i < 64; i++) {
    uint8_t * s_ = (uint8_t *) & S_log[i][0];
    for(int j=0; j < 128; j++) {
      s_[j] = rand() & 0xff;
    }
  }

  uint8_t * a_ = (uint8_t *) &A_log[0];
  for(int i=0; i<128; i++) {
    a_[i] = rand() & 0xff;
  }
}
