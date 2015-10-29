#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.c"
#include "vector.h"

v16 A[16];
v16 S_Eval[64][2][16];
//v16 Sinv_Eval[64][16];

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


#define N_BYTES 320000000

/* #define v16_cmp_eq   __builtin_ia32_pcmpeqw128 */
/* #define v16_movemask(x)   __builtin_ia32_pmovmskb128((v8) x) */

#define v16_mask(x)       (v16_movemask(x))
#define v16_nulcoef(x)    (v16_cmp_eq(x,NULL_VECT))
#define v16_pcoef(x)      (v16_cmp(x,NULL_VECT))

// Find out if a sample (v16) is invalid (i.e. if the vector has one (or more) coefficient equal to 0).

#define REJECTION_MASK(x) (v16_mask(v16_nulcoef(x)))

// Rounding of a v16 vector:

#define ROUNDING_MASK(x)  (v16_mask(v16_pcoef(x)))

#define half_bits_w(x)    (x&0x5555)
#define weak_bits8(x)     (x&0xff)
#define bits_shift_r(x,i) (x>>i)
#define bits_shift_l(x,i) (x<<i)
#define arranged_mask(x)  (bits_shift_r(half_bits_w(x), 7)^(weak_bits8(half_bits_w(x))))
#define ROUNDING(x)       (arranged_mask(ROUNDING_MASK(x)))




// res <--- a * b   (polynômes sous formes évalués)
void MultiplyPolyEval128(const v16 a[16], const v16 b[16], v16 res[16]) {
  for (int i=0; i<16; i++){
    res[i] = REDUCE_FULL_S(v16_mul(a[i], b[i]));
  }
}   

// Coef <--- coefficients du polynôme représenté sous forme évalué dans Eval
void ConvertEvalToCoefficients(const v16 Eval[16], v16 Coef[16]) {
  for(int i=0; i<16; i++) {
    Coef[i] = Eval[i];
  }

  fft128(Coef);

  // indispensable
  v16* CoefPtr = (v16*)Coef;

  for(int i=0; i<16; i++){
    CoefPtr[i] = REDUCE_FULL_S(v16_mul(CoefPtr[i], omegaPowers[i]));
  }
}


/**
 * @param x : la chaine de bits représentant la position actuelle du "gray counter"
 * @param Gray : le numéro de la valeur générée
 */
uint64_t UpdateCounterMode(uint64_t x, v16 Prod[16], const uint64_t Gray){
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
  const v16 *s = S_Eval[flip][inv];

  for(int i=0; i<16; i++) {
    Prod[i] = REDUCE_FULL_S(v16_mul(Prod[i], s[i]));
  }
  return x;
}

// renvoie le XOR de n_bytes octets de flux
unsigned char GrayCounterMode(int n_bytes){
  v16 Poly[16], Prod[16];
  uint64_t x=0, Gray_counter=0;
  int count = 0;

  unsigned char FinalOutput = 0;


  // Setup
  for(int i=0; i < 16; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes) {

    // Extraction du flux
    ConvertEvalToCoefficients(Prod, Poly);
    for(int i = 0; i < 16; i++) {
      if (REJECTION_MASK(Poly[i]) == 0) { //Perform Rejection-sampling.
        FinalOutput ^= ROUNDING(Poly[i]);
        count++;
      }
    }
    Gray_counter++;
    x = UpdateCounterMode(x, Prod, Gray_counter);
    //MultiplyPolyEval128(Eval1, Eval2, Prod); // prod' = Prod(en fait Eval1) * Eval2
  }

  return FinalOutput;
}


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
        // Sinv_Eval[i][j] = rand_v16(); // OK, Sinv n'est pas vraiment l'inverse de S. Et alors ?
      }
    }
  }
  for(int i=0; i<16; i++) {
    A[i]=rand_v16();
  }
}


int main(){
  init_secrets();

  uint64_t tsc = rdtsc();

  unsigned char Output = GrayCounterMode(N_BYTES);
  
  tsc = rdtsc() - tsc;

  printf("Output: %d\n", Output);
  printf ("%f c/B\n", 8.*tsc/(1.*N_BYTES*8));

  return 0;
}
