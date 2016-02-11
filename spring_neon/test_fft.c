#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
            
#define K 64

            
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


int test_fft128(int width, i16 omega) {
	i16 *A = malloc(N * width * sizeof(i16));
	i16 *B = malloc(N * width * sizeof(i16));

	// initialise un tableau pseudo-aléatoire
	for(int i = 0; i < N*width; i++) {
		A[i] = reduce(5*i ^ 17*i ^ 42);
	}

	// calcule les [width] FFTs parallèles en O(width * N^2).
	// B[i] = sum(A[j] * (omega^i)^j, i=0..127)
	for(int w=0; w<width; w++) {
		i16 omega_i = 1; // contient omega^i
		for(int i = 0; i < N; i++) {
			i16 omega_ij = 1; // contient omega^(ij)
			B[w + i * width] = 0;
			for(int j = 0; j < N; j++) {
				B[w + i * width] = reduce(B[w + i * width] + A[w + j*width] * omega_ij);
				omega_ij = reduce(omega_ij * omega_i);
			}
			omega_i = reduce(omega_i * omega);
		}
	}

	// check

	fft128(A);
		
	// check
	for(int i=0; i<128; i++) {
	  if (A[i] != B[i]) {
	    return 0;
	  }
	}
	return 1;

}


int main() {
  printf("parallel_reduce : %d\n", test_parallelreduce());
  printf("fft128 : %d\n", test_fft128(1,42));
  return 0; 
}
