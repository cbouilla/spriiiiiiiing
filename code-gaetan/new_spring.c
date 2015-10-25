#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

//#include "spring.h"
#include "vector.c"
#include "vector.h"
#include "poly_eval.h"
//#include "rand_generated_arrays.c"

#define CCV(x) {x, x, x, x, x, x, x, x}

const v16 NULL_VECT=CCV(0);

#define N_BYTES 69

#define v16_cmp_eq      __builtin_ia32_pcmpeqw128
#define v16_mask(x)     __builtin_ia32_pmovmskb128((v8) x)
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
#define ROUNDING(x)   (arranged_mask(ROUNDING_MASK(x)))




// multiply polynomial a and b. The result is given by res.
void MultiplyPolyEval128(v16 a[16], v16 b[16], v16 res[16]) {
  for (int i=0; i<16; i++){
    res[i] = REDUCE_FULL_S(v16_mul(a[i], b[i]));
  }
}   

void ConvertEvalToCoefficients(const v16 Eval[16], v16 Coef[16]){

  for(int i=0; i<16; i++) Coef[i]=Eval[i];

  fft128(Coef);

  v16* CoefPtr = (v16*)Coef;

  for(int i=0; i<16; i++){
    CoefPtr[i]=REDUCE_FULL_S(CoefPtr[i] * omegaPowers[i]);  // gcc transforme-t-il ça en v16_mul ?
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

flip =  __builtin_ctz(Gray);

 
 /* while((((Gray+1)>>flip)&1)==0){ */
 /*      flip++; */
 /*    } */
 /*    flip++; */

 /* if (flip != flip2) { */
 /*    printf("ca ne colle pas : gray=%" PRIx64 " / %d vs %d\n", Gray, flip, 1 + __builtin_ctz(Gray+1)); */
 /*    } */
  
 


 mask = (1 << (flip));
 x ^= mask;

 inv=(x >> flip)&1;

 for(int i=0; i<16; i++){
   Prod[i]=Prod[i]*S_Eval[inv][flip][i];
   Prod[i]=REDUCE_FULL_S(Prod[i]);
 }

 /* if(inv) { */
 /*   for(i=0; i<16; i++){ */
 /*     Prod[i]=Prod[i]*Sinv_Eval[flip][i]; */
 /*     Prod[i]=REDUCE_FULL_S(Prod[i]); */
 /*   } */
 /* } */
 /* else { */
 /*   for(i=0; i<16; i++){ */
 /*     Prod[i]=Prod[i]*S_Eval[flip][i]; */
 /*     Prod[i]=REDUCE_FULL_S(Prod[i]); */
 /*   } */
 /* } */


 return x;
  }


/** 
 * Effectue n_iterations du mode compteur en produisant le flux
 * @param n_bytes    nombre d'octets de flux à générer
 * @return           un pointeur vers une zone de mémoire contenant n_bytes octets de flux
 */
uint64_t GrayCounterMode(int n_bytes){
  v16 Poly[16], 
    //Eval2[16], 
      Prod[16];    // a * prod(s_i^x_i)
  uint64_t x=0, Gray_counter=0;
  int Output_pt=0, count=0;

  uint64_t FinalOutput=0;

  //char * Output = malloc(n_bytes);

  for(int i=0; i<16; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes) {
 

    ConvertEvalToCoefficients(Prod, Poly); //On retrouve les coefficients du polynome à partir de ses évaluations.

    for(int i=0; i<16; i++) {
      Output_pt = XOROutputUpdate(Poly[i], Output_pt, &FinalOutput);
    }
    count += Output_pt; 
   
    Gray_counter++;
    x=UpdateCounterMode(x, Prod, Gray_counter);
    //MultiplyPolyEval128(Eval1, Eval2, Prod); // prod' = Prod(en fait Eval1) * Eval2
  }
  
  return FinalOutput;
}



int main(){

  /* Test function MultiplyPolyEval128 */

  /* v16 a[16], b[16], c[16]; */

 /*  v16 A0=CCV(0), A1=CCV(1), A2=CCV(2), A3=CCV(3), A4=CCV(4), A5=CCV(5), A6=CCV(6), A7=CCV(7), A8=CCV(8), A9=CCV(9), A10=CCV(10), A11=CCV(11), A12=CCV(12), A13=CCV(13), A14=CCV(14), A15=CCV(15), B0=CCV(1), B1=CCV(1), B2=CCV(2), B3=CCV(4), B4=CCV(8), B5=CCV(16), B6=CCV(32), B7=CCV(64), B8=CCV(128), B9=CCV(256), B10=CCV(255), B11=CCV(253), B12=CCV(249), B13=CCV(241), B14=CCV(225), B15=CCV(193),  C4; */

 /*  int i; */

 /*  A0[0]=1; A1[1]=2; A2[2]=3; A3[3]=4; A4[4]=5; A5[5]=6; A6[6]=7; A7[7]=8;  */
 /*  B0[7]=7; B2[6]=6; B4[5]=5; B6[4]=4; B8[3]=3; B10[2]=2; B12[1]=1; B14[0]=0; */

 /*  a[0]=A0; a[1]=A1; a[2]=A2; a[3]=A3; a[4]=A4; a[5]=A5; a[6]=A6; a[7]=A7; a[8]=A8; a[9]=A9; a[10]=A10; a[11]=A11; a[12]=A12; a[13]=A13; a[14]=A14; a[15]=A15; */

 /*  b[0]=B0; b[1]=B1; b[2]=B2; b[3]=B3; b[4]=B4; b[5]=B5; b[6]=B6; b[7]=B7; b[8]=B8; b[9]=B9; b[10]=B10; b[11]=B11; b[12]=B12; b[13]=B13; b[14]=B14; b[15]=B15; */

  
 /*  MultiplyPolyEval128(a,b,c); */

 /* C4=c[4]; */

 /*  for(i=0; i<8; i++){ */
 /*    printf("%d x %d ----> %d\n", A4[i], B4[i], C4[i]); */
 /*  } */

  /* Test macros REJECTION_MASK and ROUNDING */
  /* v16 X; */
  /* unsigned int m1, m2, round; */

  /* X[0]=-1; X[1]=1; X[2]=-1; X[3]=1; X[4]=-1; X[5]=1; X[6]=-1; X[7]=1; */

  /* m1=REJECTION_MASK(A0); m2=REJECTION_MASK(A1); */
  /* round=ROUNDING(X); */

  /* printf("mask(A0)= %x,\n mask(A1)= %x\n round=%x\n",m1, m2, round); */

  /* Test function OutputFillIn */
 /* Test GrayCounterMode */

#ifdef rdtsc
	uint64_t tsc = rdtsc();
#else
	clock_t begin, end;
	begin = clock();
#endif

 uint64_t Output=GrayCounterMode(N_BYTES);

#ifdef rdtsc
	tsc = rdtsc() - tsc;
#else
	end = clock();
#endif


 printf("Output: %" PRIx64 "\n", Output);
  /* for(int i=0; i<30; i++) { */
  /*   printf("%02x ", (unsigned char) Output[i]); */
  /* } */

#ifdef rdtsc
	printf ("%f c/B\n", 8.*tsc/(1.*N_BYTES*8));
#else
	double dt = (double) (end - begin) / CLOCKS_PER_SEC;
	printf ("%f MB/s (time = %f)\n", ((float)N_BYTES*8/8000000dt, dt);
#endif

  printf("\n");

  return 0;
}
