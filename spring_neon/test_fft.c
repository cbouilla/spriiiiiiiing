#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
            
#define K 64
const int cste = 2;
            
#include "vector.c"

typedef short int i16;

i16 reduce(i16 a) {
	int x = (a % 257);
	if (x > 128) x -= 257;
	if (x < -128) x += 257;
	return x;
}

int test_parallelreduce() {
	v16 a, b;

	for(int i=0; i<8; i++) {
		a[i] = rand();
	}

	b = EXTRA_REDUCE(REDUCE(a));

	for(int i=0; i<8; i++) {
		if (b[i] != reduce(a[i])) {
			return 0;
		}
	}
	return 1;
}

int test_dit_butterfly() {
  v16 u, v, a, b;

  for(int i = 0; i < 8; i++) {
    u[i] = rand();
    v[i] = rand();
  }
  u = REDUCE(u);
  v = REDUCE(v);
  a = u;
  b = v;

  DIF_BUTTERFLY(u, v, cste); 
  //simulate the butterfly
  for(int i = 0; i < 8; i++){
    a[i] = a[i] + b[i];
    b[i] = a[i] - b[i] - (b[i] << cste);
  }
  for(int i = 0; i < 8; i++){
    if(a[i] !=  u[i] || b[i] != v[i]){
      return 0;
    }
  }
  return 1;
}

int main() {
  printf("parallel_reduce : %d\n", test_parallelreduce());
  printf("butterfly : %d", test_dit_butterfly());
}
