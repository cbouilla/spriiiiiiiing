#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector32.c"
#include "vector.h"


#define N_ITERATIONS 100000000

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



int main(){
  time_fft_8();
  time_fft_16();
  time_fft_128();
  return 0;
}
