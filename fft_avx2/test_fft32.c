#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
            
#define K 64            
#include "common.c"

typedef short int i16;

i16 reduce(i16 a) {
	int x = (a % 257);
	if (x > 128) x -= 257;
	if (x < -128) x += 257;
	return x;
}

int test_parallelreduce() {
	v32 a, b;

	for(int i=0; i<16; i++) {
		a[i] = rand();
	}

	b = EXTRA_REDUCE(REDUCE(a));

	for(int i=0; i<16; i++) {
		if (b[i] != reduce(a[i])) {
			return 0;
		}
	}
	return 1;
}

int test_butterfly(int k) {
	v32 a, b, c;

	for(int i = 0; i < 16; i++) {
		a[i] = rand();
	}
	a = REDUCE(a);
	b = a;

	DIT_TROUBLESOME_BUTTERFLY(a, k);

	// simulate the butterfly

	for(int i = 0; i < 8; i++) {
		c[i] = b[i] + (b[i+8] << k);
		c[i+8] = b[i] - (b[i+8] << k);
	}

	for(int i = 0; i < 16; i++) {
		if (a[i] != c[i]) {
			return 0;
		}
	}

	return 1;
}

int revbin8[] = {0, 4, 2, 6, 1, 5, 3, 7};
int revbin16[] = {0, 8, 4, 0xc, 2, 0xa, 6, 0xe, 1, 9, 5, 0xd, 3, 0xb, 7, 0xf};

int test_fft(int N, int width, i16 omega) {
	i16 *A = malloc(N * width * sizeof(i16));
	i16 *B = malloc(N * width * sizeof(i16));
	i16 *C = malloc(N * width * sizeof(i16));

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
	if (N == 8) {
		dif_fft8(A);
		// check. Output is in revbin order
		for(int i=0; i<8; i++) {
			for(int j=0; j<width; j++) {
				if (A[j + revbin8[i]*width] != B[j + i*width]) {
					return 0;
				}
			}
		}
		return 1;
	} else if (N == 16) {
		// put input in revbin order
		for(int i=0; i<16; i++)
			for(int j=0; j<8; j++)
				C[i*8 + j] = A[revbin16[i]*8 + j];

		dit_fft16(C);
		
		// check
		for(int i=0; i<16; i++) {
			for(int j=0; j<8 ;j++) {
				if (C[j + i*8] != B[j + i*8]) {
					return 0;
				}
			}
		}
		return 1;
	} else if (N == 128) {
		fft128(A);
		
		// check
		for(int i=0; i<128; i++) {
			if (A[i] != B[i]) {
				return 0;
			}
		}
		return 1;
	} else {
		printf("taille non-gérée\n");
		return 0;
	}
	return 0;
}

int test_rounding() {
	v32 A;
	uint64_t r = 0;
	// initialise un tableau pseudo-aléatoire
	for(int i = 0; i <16; i++) {
		A[i] = rand() & 0x00ff;
	}
	for(int i = 14; i >= 0; i-=2) {
		r <<= 4;
		r ^= (A[i] >> 4);

		r <<= 4;
		r ^= (A[i+1] >> 4);
	}

	uint64_t r2 = rounding(A);
	// printf("%#" PRIx64 " VS %#" PRIx64 "\n", r, r2);
	return r == r2;
}

int test_BCH() {
	v16 A;
	
	for(int i = 0; i <8; i++) {
		A[i] = rand();
	}
	
	uint64_t a = BCH128to64(A);
	uint64_t b = BCH128to64_clmul(A);
	return a == b;
}

int test_msb() {
	v32 A;
	uint16_t r = 0;
	// initialise un tableau pseudo-aléatoire
	for(int i = 0; i <16; i++) {
		A[i] = (rand() % 257) - 128;
	}
	for(int i = 15; i >= 0; i--) {
		r <<= 1;
		r ^= (A[i] > 0) ? 1 : 0;
	}

	uint16_t r2 = msb(A);
	return r == r2;
}

int test_subset_sum() {
	vLog SumA[4], SumB[4];
	init_secrets_log();
	init_subset_sum_tables();

	ComputeSubsetSum(0x123456789abcdef0, SumA);
    ComputeSubsetSum_tabulated(0x123456789abcdef0, SumB);

    for(int i=0; i<4; i++) {
    	for(int j=0; j<32; j++) {
    		if (SumA[i][j] != SumB[i][j]) {
    			return 0;
    		}
    	}
    }
    return 1;
}

int test_exponentiate_tables() {
	for(int i=0; i < 16; i++) {
		//printf("T[%d] = %d\n", i, (int16_t) ((int8_t) generatorPowersT2[i]) + 2);
	}

	for(int i=0; i<256; i++) {
		int low = i % 16;
		int high = i / 16;
		int a = (generatorPowers[i] + 257) % 257;
		int b = ((int8_t) generatorPowersT1[low] + 257) % 257;
		int c = ((int8_t) generatorPowersT2[high] + 2 + 257) % 257;

		if (a != ((b*c) % 257)) {
			return 0;
		}
	}

	return 1;
}



int test_exponentiate() {
	vLog L[4];
	v32 E1[8], E2[8];
	for(int i=0; i<4; i++) {
    	for(int j=0; j<32; j++) {
    		L[i][j] = rand() & 0xff;
    	}
	}  
   	exponentiate(L, E1);
	exponentiate_ssse3(L, E2);

	for(int a=0; a<128; a++) {
    	int i = a / 16;
    	int j = a % 16;
    	uint8_t l =  L[a / 32][a % 32];
    	if (E1[i][j] != E2[i][j]) {
	    	printf("i=%d, exp(%02x) = %d vs %d (real=%d)\n", a, l, E1[i][j], E2[i][j], generatorPowers[l]);
    		return 0;
		}
    }
    return 1;
}

int main() {
	printf("parallel_reduce : %d\n", test_parallelreduce());
	printf("butterfly : %d\n", test_butterfly(5)); 
	printf("fft8   : %d\n", test_fft(8, 16, 4)); 
	printf("fft16  : %d\n", test_fft(16, 8, 2)); 
	printf("fft128  : %d\n", test_fft(128, 1, 42)); 
	printf("rounding : %d\n", test_rounding());
	printf("BCH : %d\n", test_BCH());
	printf("msb : %d\n", test_msb());
	printf("subset_sum : %d\n", test_subset_sum());
	printf("exponentiate tables : %d\n", test_exponentiate_tables());
	printf("exponentiate : %d\n", test_exponentiate());
}
