#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#define K 64
#define EXTRACTED_BITS 4
#define N_BYTES (1024ull*1024*1024)

#include "common.c"


static const int MultiplyDeBruijnBitPosition2[32] = {
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9 };

// trick found on https://graphics.stanford.edu/~seander/bithacks.html
int log_for_powers_of_two(const uint32_t x) {
  return MultiplyDeBruijnBitPosition2[(uint32_t)(x * 0x077CB531u) >> 27];
}
 

// renvoie le XOR de n_bytes octets de flux
uint64_t OutputFeedbackMode(uint64_t n_bytes){
  vLog Sum[4];
  v32 Poly[8];
  uint64_t O[8];
  uint64_t count = 0;
  uint64_t best = 0, good=0, bad=0;
  uint64_t FinalOutput = 0;
  int ptr = 0;     // the next available bits must be appended to O[ptr]
  int spaces = 64; // how many bits must be appended to O[ptr]

  O[0] = 0;
  O[1] = 0;

  while(count < n_bytes) {
    
    //printf("%016" PRIx64 " \n", FinalOutput);
    ComputeSubsetSum_tabulated(O[0], Sum);

    // Extraction du flux
    ConvertSubsetSumToCoefficients(Sum, Poly);
    
    ptr = 0;
    spaces = 64;
    for(int i = 0; i < 8; i++) {
      O[i] = 0;
    }

    for(int i = 0; i < 8; i++) {
      const v32 a = Poly[i];
      const uint64_t x = rounding(a); // all nibbles may not necessarily be OK
      const uint32_t r = reject(a);
      if ((spaces == 64) && (r == 0)) { // easy case: 64 bits fit into 64 bits. This is a (faster) special case of the complicated code below.
        O[ptr] = x;
        ptr++;
        best++;
      } else { // annoying case. Either we are missing a nibble now, or we missed one/a few before.

        uint64_t remaining = x;  // holds the reconstructed value
        int remaining_size = 64; // size in bits

        if (r != 0) { // missing nibble(s). Eliminate them.
          uint32_t mask = r & 0x55555555;

          if ((mask & (mask - 1)) == 0) {  // just one nibble is missing (most frequent case)
            int k = log_for_powers_of_two(mask) * 2; 
            uint64_t selector = (1 << k) - 1;
            remaining = (x & selector) ^ ((x >> 4) & (~selector));
            remaining_size = 60;
            good++;
          } else { // more than one nibble is missing. Slow, generic procedure
            remaining = 0;
            remaining_size = 0;
            for(int k=0; k<16; k++) { // we check each nibble individually
              if ((mask & 1) == 0) { // bottom nibble is OK
                uint64_t nibble = (x >> (4ull * k)) & 0xf;
                remaining ^= nibble << remaining_size;
                remaining_size += 4;
              }
              mask >>= 2;
            }
            bad++;
          } 
        } // end building "remaining"


        // append "remaining" to the output
        if (remaining_size > spaces) {  // remaining does not fit into O[ptr]
          O[ptr] ^= (remaining & ((1ull << spaces) - 1)) << (64ull - spaces);
          remaining >>= spaces;
          remaining_size -= spaces;
          spaces = 64;
          ptr++;
        }

        // now it fits entirely
        O[ptr] ^= remaining << (64ull - spaces);
        spaces -= remaining_size;
        if (spaces == 0) {
          spaces = 64;
          ptr++;
        }
      } // end annoying case
      

      if (ptr >= 6) { // exit the loop as soon as enough output is ready
        break;
      }
    } //end for
    
    FinalOutput ^= O[0];
    FinalOutput ^= O[1];
    FinalOutput ^= O[2];
    FinalOutput ^= O[3];
    FinalOutput ^= O[4];
    FinalOutput ^= O[5];

    count += 48;
  }
  printf("best=%llu, good=%llu, bad=%llu\n", best, good, bad);
  return FinalOutput;
}


int main(){
  init_secrets_log();
  init_subset_sum_tables();

  uint64_t tsc = rdtsc();

  uint64_t Output = OutputFeedbackMode(N_BYTES);
  
  tsc = rdtsc() - tsc;

  printf("Output: %#" PRIx64 "\n", Output);
  printf("%f c/B\n", tsc/(1.*N_BYTES));

  return 0;
}
