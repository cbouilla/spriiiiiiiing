#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "common.c"

#define EXTRACTED_BITS 4
#define N_BYTES 320000000

// renvoie le XOR de n_bytes octets de flux
uint64_t OutputFeedbackMode(int n_bytes){
  vLog Sum[4];
  v32 Poly[8];
  uint64_t x = 0;
  int count = 0;

  uint64_t FinalOutput = 0;


  while(count < n_bytes) {
    ComputeSubsetSum_tabulated(x, Sum);

    // Extraction du flux
    ConvertSubsetSumToCoefficients(Sum, Poly);
    uint16_t M[8];
    for(int i = 0; i < 8; i++) {
      M[i] = msb(Poly[i]);
    }
    // apply the BCH code
    v16 biased = {M[0], M[1], M[2], M[3], M[4], M[5], M[6], M[7]};
    x = BCH128to64_clmul(biased);
    FinalOutput ^= x;

    count += 8;
  }

  return FinalOutput;
}


int main(){
  init_secrets_log();
  init_subset_sum_tables();

  uint64_t tsc = rdtsc();

  uint64_t Output = OutputFeedbackMode(N_BYTES);
  
  tsc = rdtsc() - tsc;

  printf("Output: %#" PRIx64 "\n", Output);
  printf ("%f c/B\n", tsc/(1.*N_BYTES));

  return 0;
}
