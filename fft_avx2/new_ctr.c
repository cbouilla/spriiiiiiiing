#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#define K 64
#define N_BYTES (1024ull*1024*1024)

#include "common.c"

// renvoie le XOR de n_bytes octets de flux
uint64_t GrayCounterMode(int n_bytes){
  v32 Poly[8], Prod[8];
  uint64_t x=0, Gray_counter=0;
  int count = 0;

  uint64_t FinalOutput = 0;

  // Setup
  for(int i=0; i < 8; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes) {
    ConvertEvalToCoefficients(Prod, Poly);
    for(int i = 0; i < 8; i++) {
      const v32 a = Poly[i];
      if (!reject(a)) {
          FinalOutput ^= rounding(a);
          count += 8;
      }
    }
    Gray_counter++;
    x = UpdateCounterMode(x, Prod, Gray_counter);
  }

  return FinalOutput;
}


int main(){
  init_secrets();

  uint64_t tsc = rdtsc();

  uint64_t Output = GrayCounterMode(N_BYTES);
  
  tsc = rdtsc() - tsc;

  printf("Output: %#" PRIx64 "\n", Output);
  printf ("%f c/B\n", tsc/(1.*N_BYTES));

  return 0;
}
