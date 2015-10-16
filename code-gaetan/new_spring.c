#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

//#include "spring.h"
#include "vector.c"
#include "vector.h"
//#include "rand_generated_arrays.c"

#define CCV(x) {x, x, x, x, x, x, x, x}

const v16 NULL_VECT=CCV(0);

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
  v16 A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15;

  A0=a[0];  A1=a[1];  A2=a[2];  A3=a[3];  A4=a[4];  A5=a[5];  A6=a[6];  A7=a[7]; A8=a[8];  A9=a[9];  A10=a[10];  A11=a[11];  A12=a[12];  A13=a[13];  A14=a[14];  A15=a[15];
 
B0=b[0];  B1=b[1];  B2=b[2];  B3=b[3];  B4=b[4];  B5=b[5];  B6=b[6];  B7=b[7];  B8=b[8];  B9=b[9];  B10=b[10];  B11=b[11];  B12=b[12];  B13=b[13];  B14=b[14];  B15=b[15];

 X0 = v16_mul(A0,B0);
 X1 = v16_mul(A1,B1);
 X2 = v16_mul(A2,B2);
 X3 = v16_mul(A3,B3);
 X4 = v16_mul(A4,B4);
 X5 = v16_mul(A5,B5);
 X6 = v16_mul(A6,B6);
 X7 = v16_mul(A7,B7);
 X8 = v16_mul(A8,B8);
 X9 = v16_mul(A9,B9);
 X10 = v16_mul(A10,B10);
 X11 = v16_mul(A11,B11);
 X12 = v16_mul(A12,B12);
 X13 = v16_mul(A13,B13);
 X14 = v16_mul(A14,B14);
 X15 = v16_mul(A15,B15);

#define X(i) X##i

DO_REDUCE_FULL_S(0);
DO_REDUCE_FULL_S(1);
DO_REDUCE_FULL_S(2);
DO_REDUCE_FULL_S(3);
DO_REDUCE_FULL_S(4);
DO_REDUCE_FULL_S(5);
DO_REDUCE_FULL_S(6);
DO_REDUCE_FULL_S(7);
DO_REDUCE_FULL_S(8);
DO_REDUCE_FULL_S(9);
DO_REDUCE_FULL_S(10);
DO_REDUCE_FULL_S(11);
DO_REDUCE_FULL_S(12);
DO_REDUCE_FULL_S(13);
DO_REDUCE_FULL_S(14);
DO_REDUCE_FULL_S(15);

 res[0]=X0;  res[1]=X1;  res[2]=X2;  res[3]=X3;  res[4]=X4;  res[5]=X5;  res[6]=X6;  res[7]=X7;  res[8]=X8;  res[9]=X9;  res[10]=X10;  res[11]=X11;  res[12]=X12;  res[13]=X13;  res[14]=X14;  res[15]=X15;

}   


uint64_t OutputFillIn(v16 a, int* i, uint64_t Output){
  
  //Perform Rejection-sampling.
  if(REJECTION_MASK(a)) return Output;

  uint64_t tmp;

  //Round vector.
  tmp=bits_shift_l(ROUNDING(a), 8*(*i));
  Output=Output^tmp;
  (*i)++;

  return Output;
}


int main(){

  /* Test function MultiplyPolyEval128 */

  v16 a[16], b[16], c[16];

  v16 A0=CCV(0), A1=CCV(1), A2=CCV(2), A3=CCV(3), A4=CCV(4), A5=CCV(5), A6=CCV(6), A7=CCV(7), A8=CCV(8), A9=CCV(9), A10=CCV(10), A11=CCV(11), A12=CCV(12), A13=CCV(13), A14=CCV(14), A15=CCV(15), B0=CCV(1), B1=CCV(1), B2=CCV(2), B3=CCV(4), B4=CCV(8), B5=CCV(16), B6=CCV(32), B7=CCV(64), B8=CCV(128), B9=CCV(256), B10=CCV(255), B11=CCV(253), B12=CCV(249), B13=CCV(241), B14=CCV(225), B15=CCV(193), C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15;

  int i;

  A0[0]=1; A1[1]=2; A2[2]=3; A3[3]=4; A4[4]=5; A5[5]=6; A6[6]=7; A7[7]=8; 
  B0[7]=7; B2[6]=6; B4[5]=5; B6[4]=4; B8[3]=3; B10[2]=2; B12[1]=1; B14[0]=0;

  a[0]=A0; a[1]=A1; a[2]=A2; a[3]=A3; a[4]=A4; a[5]=A5; a[6]=A6; a[7]=A7; a[8]=A8; a[9]=A9; a[10]=A10; a[11]=A11; a[12]=A12; a[13]=A13; a[14]=A14; a[15]=A15;

  b[0]=B0; b[1]=B1; b[2]=B2; b[3]=B3; b[4]=B4; b[5]=B5; b[6]=B6; b[7]=B7; b[8]=B8; b[9]=B9; b[10]=B10; b[11]=B11; b[12]=B12; b[13]=B13; b[14]=B14; b[15]=B15;

  
  MultiplyPolyEval128(a,b,c);

  C0=c[0]; C1=c[1]; C2=c[2]; C3=c[3]; C4=c[4]; C5=c[5]; C6=c[6]; C7=c[7]; C8=c[8]; C9=c[9]; C10=c[10]; C11=c[11]; C12=c[12]; C13=c[13]; C14=c[14]; C15=c[15];

  for(i=0; i<8; i++){
    printf("%d x %d ----> %d\n", A4[i], B4[i], C4[i]);
  }

  /* Test macros REJECTION_MASK and ROUNDING */
  v16 X;
  unsigned int m1, m2, round;

  X[0]=-1; X[1]=1; X[2]=-1; X[3]=1; X[4]=-1; X[5]=1; X[6]=-1; X[7]=1;

  m1=REJECTION_MASK(A0); m2=REJECTION_MASK(A1);
  round=ROUNDING(X);

  printf("mask(A0)= %x,\n mask(A1)= %x\n round=%x\n",m1, m2, round);

  /* Test function OutputFillIn */
  uint64_t Output=0;
  i=0;

  Output=OutputFillIn(X, &i, Output);
  printf("Output : %x \n", (unsigned) Output);
  Output=OutputFillIn(A0, &i, Output);
  printf("Output : %x \n", (unsigned) Output);
  Output=OutputFillIn(X, &i, Output);
  printf("Output : %x \n", (unsigned) Output);
  printf("%d\n", i);

  return 0;
}
