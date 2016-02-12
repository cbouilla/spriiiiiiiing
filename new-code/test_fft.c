#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
            
#define K 64

#include "vector.h"           
#include "vector.c"

typedef short int i16;

i16 reduce(i16 a) {
	int x = (a % 257);
	if (x > 128) x -= 257;
	if (x < -128) x += 257;
	return x;
}


int test_fft_64(int N, int width, i16 omega) {
 i16 A[N*width] __attribute ((aligned (16)));
  i16 B[N*width];
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

        
	fft64(A);

	//check
	for(int i = 0; i < 64; i++){
	  printf("A[%d] = %d B = %d\n", i, A[i], B[i]);
	}

	for(int i = 0; i < 64; i++){
	  if(A[i] != B[i]) {
	    printf("A[%d] = %d vs B[%d] = %d\n", i, A[i], i, B[i]);

	    return 0;
	  }
	}
	return 1;


}


int main(){
printf("fft64 : %d\n", test_fft_64(64, 1, -35));
}
