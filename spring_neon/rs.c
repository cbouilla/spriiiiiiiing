#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#define K 64
#define N_BYTES (1024ull*1024*1024)

#include "common.c"

uint32_t GrayCounterMode(int n_bytes){
  v16 Poly[16], Prod[16];
  uint32_t x=0, Gray_Counter=0;
  int count = 0;

  uint32_t FinalOutput = 0;

  //Setup :

  for(int i = 0; i < 16; i++){
    Prod[i] = A[i];
  }

  while(count < n_bytes) {
    ConvertEvalToCoefficients(Prod, Poly);
    for(int i = 0; i < 8; i++){
      const v16 a = Poly[i];
      if(!reject(a)) {
	FinalOutput ^= rounding4(a);
	count += NBITS;
      }
    }
    Gray_Counter++;
    x = UpdateCounterMode(x, Prod, Gray_Counter);
  }
  return FinalOutput;
}


int main(){

  clock_t begin, end;
  begin = clock();


  init_secrets();  // TODO remplacer par un LFSR ?

  begin = clock();

  uint32_t Output = GrayCounterMode(N_BYTES);

  end = clock();


  printf("Output : %#" PRIx32 "\n", Output);


  double dt = (double) (end - begin) / CLOCKS_PER_SEC;
  printf ("%f MB/s (time = %f)\n", ((float)N_BYTES/1000000)/dt, dt);

  return 0;
}
