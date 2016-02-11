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
	i16 *A = malloc(128 * width * sizeof(i16));
	i16 *B = malloc(128 * width * sizeof(i16));

	printf("DEBUG : L42\n");

	// initialise un tableau pseudo-aléatoire
	for(int i = 0; i < 128*width; i++) {
		A[i] = reduce(5*i ^ 17*i ^ 42);
	}

	printf("DEBUG : L49\n");

	// calcule les [width] FFTs parallèles en O(width * N^2).
	// B[i] = sum(A[j] * (omega^i)^j, i=0..127)
	for(int w=0; w<width; w++) {
		i16 omega_i = 1; // contient omega^i
		for(int i = 0; i < 128; i++) {
			i16 omega_ij = 1; // contient omega^(ij)
			B[w + i * width] = 0;
			printf("DEBUG : L58; w %d, i %d\n", w, i);

			for(int j = 0; j < 128; j++) {
				B[w + i * width] = reduce(B[w + i * width] + A[w + j*width] * omega_ij);

				printf("DEBUG : L63; j %d\n", j);

				omega_ij = reduce(omega_ij * omega_i);
			}
			omega_i = reduce(omega_i * omega);
		}
	}

	// check
	printf("DEBUG : L72\n");
	fft128(A);
	printf("DEBUG : L74\n");
		
	// check
	for(int i=0; i<128; i++) {
	  printf("DEBUG : L78 i %d\n", i);
	  if (A[i] != B[i]) {
	    return 0;
	  }
	  printf("DEBUG : L82\n");
	}
	return 1;

}


int main() {
  printf("parallel_reduce : %d\n", test_parallelreduce());
  printf("fft128 : %d\n", test_fft128(1,42));
  return 0; 
}
