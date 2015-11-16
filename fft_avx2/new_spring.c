#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector32.c"
#include "vector.h"

#define EXTRACTED_BITS 4
#define N_BYTES 320000000


v32 A[8];
v32 S_Eval[64][2][8];
//v16 Sinv_Eval[64][16];

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


const v32 NULL_VECT = v32_cst(0); 
const v32 BIT6 = v32_cst(4095);
const v32 BIT5 = v32_cst(1024);
const v32 BIT4 = v32_cst(255);





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

// renvoie le XOR de n_bytes octets de flux
unsigned short int GrayCounterMode(int n_bytes){
  v32 Poly[8], Prod[8];
  uint64_t x=0, Gray_counter=0;
  int count = 0;

  unsigned short FinalOutput = 0;


  // Setup
  for(int i=0; i < 8; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes) {

    // Extraction du flux
    ConvertEvalToCoefficients(Prod, Poly);
    for(int i = 0; i < 8; i++) {
      v32 Y;
      if (v32_movemask(v32_cmp_eq(Poly[i], NULL_VECT)) == 0) { //Perform Rejection-sampling.
          /*if (EXTRACTED_BITS > 1) {
            Y = Poly[i] * Poly[i];
          }*/

          // extract signs
          int m = v32_movemask(v32_cmp_gt(Poly[i], NULL_VECT));

          if (EXTRACTED_BITS == 1) {
            FinalOutput ^=  ((m & 0x55550000) >> 15) ^ (m & 0x00005555);
            count += 2;
            continue;
          }

          // square the vector to get absolute values
          Y = Poly[i] * Poly[i];
          int n = v32_movemask(v32_cmp_gt(Y, BIT6));
          short int high = (n & 0x55550000) ^ (m & 0xcccc0000);
          short int low = (n & 0x00005555) ^ (m & 0x0000cccc);

          FinalOutput ^= high >> 16;
          FinalOutput ^=  low;

          if (EXTRACTED_BITS == 2) {  
            count += 4;
            continue;
          }
  
          int o = v32_movemask(v32_cmp_gt(Y, BIT5));
          if (EXTRACTED_BITS == 3) {
            FinalOutput ^=  ((o & 0x55550000) >> 15) ^ (o & 0x00005555);
            count += 6;
            continue;
          }

          int p = v32_movemask(v32_cmp_gt(Y, BIT4));
          short int third =  (o & 0x55550000) ^ (p & 0xcccc0000);
          short int fourth = (o & 0x00005555) ^ (p & 0x0000cccc);

          FinalOutput ^= third >> 16;
          FinalOutput ^= fourth;
          count += 8;
      }
    }
    Gray_counter++;
    x = UpdateCounterMode(x, Prod, Gray_counter);
  }

  return FinalOutput;
}


v32 rand_v32() {
  v32 x;
  for(int j=0; j < 16; j++) {
    do {
      x[j] = rand();
    } while(x[j] == 0);
  }
  return REDUCE_FULL(x);
}

// méthode top-secrète pour initialiser A et les s_i. Ne pas divulguer au public !
void init_secrets() {
  srand(42);

  for(int i=0; i < 64; i++) {
    for(int k=0; k<2; k++){
      for(int j=0; j < 8; j++) {
        S_Eval[i][k][j] = rand_v32();
      }
    }
  }
  for(int i=0; i<8; i++) {
    A[i] = rand_v32();
  }
}


int main(){
  init_secrets();

  uint64_t tsc = rdtsc();

  unsigned short int Output = GrayCounterMode(N_BYTES);
  
  tsc = rdtsc() - tsc;

  printf("Output: %d\n", Output);
  printf ("%f c/B\n", 8.*tsc/(1.*N_BYTES*8));

  return 0;
}
