#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>


#include "vector.c"
#include "vector.h"

#define N_ITERATIONS 1000000

extern uint64_t BCH128to64(uv64 in);

uint64_t rdtsc() {
  return clock() * 900ull;
}

void time_fft8() {
   v16 stuff[8];

  // deterministic initialization
  srand(42);
  for(int i=0; i<8; i++) {
    for(int j=0; j < 8; j++) {
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
    stuff[i] = REDUCE_FULL(stuff[i]);
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
    for(int j=0; j < 8; j++)
      printf("%04d\t", stuff[i][j]);
    printf("\n");
  }

  printf ("1x FFT-128 : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
  printf ("= %f cycles/byte\n", tsc/(1.*N_ITERATIONS)/64);
}

void time_BCH() {
  uv64 x = {0x12345678 + rdtsc(), 0xabcdef + rdtsc()};
  
  // warm-up
  uint64_t y = BCH128to64(x);

  uint64_t tsc = rdtsc();
  for(int i=0; i< N_ITERATIONS; i++) {
    y += BCH128to64(x);
  } 
  tsc = rdtsc() - tsc;

  printf("%lld\n", y);

  printf ("BCH code : %f cycles/iterations\n", tsc/(1.*N_ITERATIONS));
  printf ("= %f cycles/byte\n", tsc/(1.*N_ITERATIONS)/8);


}

int main(){
  //time_fft8();
  //time_fft64();
  //time_fft128();
  time_BCH();
  return 0;
}
