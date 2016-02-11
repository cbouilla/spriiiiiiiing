#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.h"

v16 A[16];
v16 S_Eval[64][2][16];

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

const v16 NULL_VECT = {0, 0, 0, 0, 0, 0, 0, 0};
const v16 VECT128 = {128, 128, 128, 128, 128, 128, 128, 128};

#define N_BYTES 320000000
#define N_BITS 4

v16 rand_v16() {
  v16 x;
  for(int j=0; j < 16; j++) {
    do {
      x[j] = rand();
    } while(x[j]==0);
  }
  return REDUCE_FULL_S(x);
}

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
