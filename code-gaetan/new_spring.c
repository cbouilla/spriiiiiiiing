#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.c"
#include "vector.h"
#include "poly_eval.h"


#define CCV(x) {x, x, x, x, x, x, x, x}

const v16 NULL_VECT = CCV(0);

#define N_BYTES 80000000

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
 * @param a : le vecteur à traiter
 * @param i : la position du prochain octet à stocker dans Output
 * @param Ouutput : le flux de sortie du PRNG
 * @return position (éventuellement modifiée) du prochain octet à stocker dans Output
 */
int OutputFillIn(v16 a, int i, char *Output){
  
  //Perform Rejection-sampling.
  if(REJECTION_MASK(a)) {
    return i;
  }
  Output[i] = ROUNDING(a);
  return i + 1;
}


/**
 * @param a : vecteur à traiter
 * @param i : Position du prochain octet, mod 8.
 * @param Output : Sortie du PRNG Xorée avec elle-meme.
 * @return position du prochain octet à stocker mod 8.
 */
int XOROutputUpdate(v16 a, int i, uint64_t *Output){
  assert (i<8);

    if(REJECTION_MASK(a)) {
      return i;
    }

  uint64_t tmp;
  tmp=ROUNDING(a);
  tmp=tmp<<(8*i);
  (*Output) ^= tmp;
  i++;

  return (i&0x7);
}

/**
 * @param x : la chaine de bits représentant la position actuelle du "gray counter"
 * @param Gray : le numéro de la valeur générée
 */
uint64_t UpdateCounterMode(uint64_t x, v16 Prod[16], const uint64_t Gray){
  int flip;
  uint64_t mask, inv;

  // index du bit qui doit changer
  flip = 1 + __builtin_ctz(Gray+1);

  mask = (1 << (flip));
  x ^= mask;

  // nouvelle valeur du bit qui a changé
  inv = x & mask;

  // jeu de s_i à utiliser
  const v16 *s = inv ? Sinv_Eval[flip] : S_Eval[flip];

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
  for(int i=0; i<16; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes) {

    // Extraction du flux
    ConvertEvalToCoefficients(Prod, Poly);
    for(int i=0; i<16; i++) {
      if (REJECTION_MASK(Poly[i]) == 0) { //Perform Rejection-sampling.
        FinalOutput ^= ROUNDING(Poly[i]);
        count++;
      }
    }

    // mise à jour de l'état interne
    x = UpdateCounterMode(x, Prod, Gray_counter);
    Gray_counter++;
  }

  return FinalOutput;
}


v16 rand_v16() {
  v16 x;
  for(int j=0; j < 16; j++) {
    x[j] = rand();
  }
  return REDUCE_FULL_S(x);
}

// méthode top-secrète pour initialiser A et les s_i. Ne pas divulguer au public !
void init_secrets() {
  srand(42);

  for(int i=0; i < 64; i++) {
    for(int j=0; j < 16; j++) {
      S_Eval[i][j] = rand_v16();
      Sinv_Eval[i][j] = rand_v16(); // OK, Sinv n'est pas vraiment l'inverse de S. Et alors ?
    }
    A[i] = rand_v16();
  }
}


int main(){
  init_secrets();

#ifdef rdtsc
	uint64_t tsc = rdtsc();
#else
	clock_t begin, end;
	begin = clock();
#endif

  unsigned char Output = GrayCounterMode(N_BYTES);

#ifdef rdtsc
	tsc = rdtsc() - tsc;
#else
	end = clock();
#endif

 printf("Output: %d\n", Output);
  
#ifdef rdtsc
	printf ("%f c/B\n", 8.*tsc/(1.*N_BYTES*8));
#else
	double dt = (double) (end - begin) / CLOCKS_PER_SEC;
	printf ("%f MB/s (time = %f)\n", ((float)N_BYTES*8/8000000dt, dt);
#endif

  return 0;
}
