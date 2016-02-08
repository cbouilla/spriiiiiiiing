#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#define K 64
#include "common.c"


#define N_ITERATIONS 100000001

void time_fft_8() {
  v32 stuff[8];

  // deterministic initialization
  srand(42);
  for(int i=0; i<8; i++) {
    for(int j=0; j < 16; j++) {
      stuff[i][j] = rand();
    }
    stuff[i] = REDUCE_FULL(stuff[i]);
  }

  // warm-up
  dif_fft8(stuff);

  uint64_t tsc = rdtsc();

  for(int i=0; i < N_ITERATIONS; i++) {
    dif_fft8(stuff);
  }

  tsc = rdtsc() - tsc;

  // display to force compiler not to skip the fft (+ human checking...)
  for(int i=0; i<8; i++) {
    printf("x[%02d] = ", i);
    for(int j=0; j < 16; j++)
      printf("%04d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("16x parallel DIF FFT-8 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}


void time_fft_16() {
  v32 stuff[8];

  // deterministic initialization
  srand(42);
  for(int i=0; i<8; i++) {
    for(int j=0; j < 16; j++) {
      stuff[i][j] = rand();
    }
    stuff[i] = REDUCE_FULL(stuff[i]);
  }

  // warm-up
  dit_fft16(stuff);

  uint64_t tsc = rdtsc();

  for(int i=0; i < N_ITERATIONS; i++) {
    dit_fft16(stuff);
  }

  tsc = rdtsc() - tsc;

  // display to force compiler not to skip the fft (+ human checking...)
  for(int i=0; i<8; i++) {
    printf("x[%02d] = ", i);
    for(int j=0; j < 16; j++)
      printf("%04d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("8x parallel DIT FFT-16 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}

void time_fft_128() {
  v32 stuff[8];

  // deterministic initialization
  srand(42);
  for(int i=0; i<8; i++) {
    for(int j=0; j < 16; j++) {
      stuff[i][j] = rand();
    }
    stuff[i] = REDUCE_FULL(stuff[i]);
  }

  // warm-up
  fft128(stuff);

  uint64_t tsc = rdtsc();

  for(int i=0; i < N_ITERATIONS; i++) {
    fft128(stuff);
  }

  tsc = rdtsc() - tsc;

  // display to force compiler not to skip the fft (+ human checking...)
  for(int i=0; i<8; i++) {
    printf("x[%02d] = ", i);
    for(int j=0; j < 16; j++)
      printf("%4d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("1x FFT-128 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}


void time_subset_sum() {
  vLog Sum[4], Total[4];

  Total[0] = ZERO_VECT;
  Total[1] = ZERO_VECT;
  Total[2] = ZERO_VECT;
  Total[3] = ZERO_VECT;

  uint64_t tsc = rdtsc();
  for(int i=0; i < N_ITERATIONS; i++) {
    ComputeSubsetSum(0x123456789abcdef0 + i, Sum);
    //ComputeSubsetSum_tabulated(0x123456789abcdef0 + i, Sum);
    Total[0] += Sum[0];
    Total[1] += Sum[1];
    Total[2] += Sum[2];
    Total[3] += Sum[3];
  }
  tsc = rdtsc() - tsc;

  // display to force compiler not to skip the fft (+ human checking...)
  for(int i=0; i<4; i++) {
    printf("x[%02d] = ", i);
    for(int j=0; j < 32; j++)
      printf("%02x ", (unsigned char) Total[i][j]);
    printf("\n");
  }

  printf ("1x subset-sum  : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}
  
void time_exponentiate() {
  vLog Sum[4];
  v32 Total[8], Coef[8];

  Sum[0] = A_log[0];
  Sum[1] = A_log[1];
  Sum[1] = A_log[2];
  Sum[2] = A_log[3];

  for(int i=0; i<8; i++) {
    Total[i] = ZERO_VECT;
  }
  
  uint64_t tsc = rdtsc();
  for(int i=0; i < N_ITERATIONS; i++) {
    exponentiate(Sum, Coef);
    //ComputeSubsetSum_tabulated(0x123456789abcdef0 + i, Sum);
  
    for(int i=0; i<8; i++) {
      Total[i] += Coef[i];
    }
  
  }
  tsc = rdtsc() - tsc;

  // display to force compiler not to skip the fft (+ human checking...)
  for(int i=0; i<4; i++) {
    printf("x[%02d] = ", i);
    for(int j=0; j < 32; j++)
      printf("%02x ", (unsigned char) Total[i][j]);
    printf("\n");
  }

  printf ("1x subset-sum  : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}

int main(){
  init_secrets_log();
 /* time_fft_8();
  time_fft_16();
  time_fft_128();
  init_subset_sum_tables()
  time_subset_sum(); */
  time_exponentiate();

  return 0;
}
