#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>

#define K 64
#define N_BYTES (1024ull*1024*1024)

#include "common.c"

// renvoie le XOR de n_bytes octets de flux
uint64_t GrayCounterMode(int n_bytes){
  v16 Poly[16], Prod[16];

  uint32_t x=0, Gray_counter=0;
  int count = 0;

  uint32_t FinalOutput = 0;
  uint64_t output;

  // Setup
  for(int i=0; i < 16; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes){
    //Extraction du flux
    ConvertEvalToCoefficients(Prod, Poly);
    char M[16];
    for(int i =0; i < 16; i++){
      M[i] = PermutedMSB(Poly[i]);
    }
    // apply BCH code
    v8 biased = {M[0], M[1], M[2], M[3], M[4], M[5], M[6], M[7], M[8], M[9], M[10], M[11], M[12], M[13], M[14], M[15]};

    output = BCH128to64(vreinterpretq_u64_s8(biased));
    FinalOutput ^= (output & 0x00000000ffffffff);
    FinalOutput ^= (output >> 32);
    count +=8;
    Gray_counter++;
    x = UpdateCounterMode(x, Prod, Gray_counter);
  }

  return FinalOutput;
}

int main(){
  init_secrets();

  clock_t begin, end;
  begin = clock();

  uint32_t Output = GrayCounterMode(N_BYTES);

  end = clock();

  printf("Output :%"PRIx32"\n", Output);

  double dt = (double) (end - begin) / CLOCKS_PER_SEC;
  printf ("%f MB/s (time = %f)\n", ((float)N_BYTES/1000000)/dt, dt);

  printf ("%.2f cycles / byte\n", dt*1e9 / N_BYTES);

  return 0;
}
