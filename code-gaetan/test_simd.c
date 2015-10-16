
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
            
#include "vector.h"
#include "vector.c"

typedef short int i16;

i16 reduce(i16 a) {
	int x = (a % 257);
	if (x > 128) x -= 257;
	if (x < -128) x += 257;
	return x;
}


int test_inverse_fft (int N, i16 omega, i16 final_coeff) {
	i16 *A = malloc(N * sizeof(i16));
	i16 *B = malloc(N * sizeof(i16));

	// initialise un tableau pseudo-aléatoire
	srand(42);
	for(int i = 0; i < N; i++) {
		A[i] = reduce(rand());
	}

	// calcule la FFT en O(n^2).
	// B[i] = sum(A[j] * (omega^i)^j, i=0..127)
	i16 omega_i = 1; // contient omega^i
	for(int i = 0; i < N; i++) {
		i16 omega_ij = 1; // contient omega^(ij)
		B[i] = 0;
		for(int j = 0; j < N; j++) {
			B[i] = reduce(B[i] + A[j] * omega_ij);
			omega_ij = reduce(omega_ij * omega_i);
		}
		omega_i = reduce(omega_i * omega);
	}

	// maintenant, calcule l'opération réciproque...

	if (N == 64) {
		fft64(B);
	} else if (N == 128) {
		fft128(B);
	} else {
		printf("taille non-gérée\n");
		return 0;
	}

	// (opération supplémentaire : multiplication par 1/N)
	for(int i=0; i < N; i++) {
		B[i] = reduce(B[i] * final_coeff);
	}

	// ...puis vérifie qu'on retombre bien sur le truc de départ
	for(int i=0; i < N; i++) {
		if (A[i] != B[i]) {
			return 0;
		}
	}
	return 1;
}


int main() {
	assert (test_inverse_fft(64, 95, -4) == 1); // 46^(-1) = 94, 64^(-1) = 125
	assert (test_inverse_fft(128, 98, -2) == 1);
}