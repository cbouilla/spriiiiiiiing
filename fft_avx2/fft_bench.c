#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "vector.c"
#include "vector.h"


#define N_ITERATIONS 100000000

void time_fft8() {
   v16 stuff[8];

  // deterministic initialization
  srand(42);
  for(int i=0; i<8; i++) {
    for(int j=0; j < 8; j++) {
      stuff[i][j] = rand();
    }
    stuff[i] = REDUCE_FULL_S(stuff[i]);
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
    for(int j=0; j < 8; j++)
      printf("%04d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("8x parallel DIF FFT-8 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}

void time_fft64() {
   v16 stuff[8];

  // deterministic initialization
  srand(42);
  for(int i=0; i<8; i++) {
    for(int j=0; j < 8; j++) {
      stuff[i][j] = rand();
    }
    stuff[i] = REDUCE_FULL_S(stuff[i]);
  }

  // warm-up
  fft64(stuff);

  uint64_t tsc = rdtsc();

  for(int i=0; i < N_ITERATIONS; i++) {
    fft64(stuff);
  }

  tsc = rdtsc() - tsc;

  // display to force compiler not to skip the fft (+ human checking...)
  for(int i=0; i<8; i++) {
    printf("x[%02d] = ", i);
    for(int j=0; j < 8; j++)
      printf("%04d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("1x FFT-64 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}

void time_fft128() {
   v16 stuff[16];

  // deterministic initialization
  srand(42);
  for(int i=0; i<16; i++) {
    for(int j=0; j < 8; j++) {
      stuff[i][j] = rand();
    }
    stuff[i] = REDUCE_FULL_S(stuff[i]);
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
    for(int j=0; j < 8; j++)
      printf("%04d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("1x FFT-128 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
}

int main(){
  time_fft8();
  time_fft64();
  time_fft128();
  return 0;
}
