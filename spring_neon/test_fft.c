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


int test_dif_butterfly() {
  v16 a, b, c, d;
  int tmp;

  for(int i = 0; i < 8; i++){
    a[i] = rand();
    b[i] = rand();
  }
  a = REDUCE(a);
  b = REDUCE(b);
  c = a;
  d = b;

 
  DIF_BUTTERFLY(a, b, 1);

 
  //simulate the butterfly

  for(int i = 0; i < 8; i++) {
    tmp = c[i];
    c[i] = c[i] + d[i];
    d[i] = ((tmp - d[i]) << 2);
  }

  for(int i = 0; i < 8; i++){
    if(a[i] != c[i] ){
      printf("i : %d, a : %d, c : %d\n", i, a[i], c[i]);
      return 0;
    }
  }

  for(int i = 0; i < 8; i++){
    if(b[i] != d[i]){
      printf("i : %d, b : %d, d : %d\n", i, b[i], d[i]);
      return 0;
    }
  }

  return 1;
}

int revbin8[] = {0, 4, 2, 6, 1, 5, 3, 7};

int test_fft8(int width, i16 omega) {
  i16 A[8*width] __attribute ((aligned (16)));
  i16 B[8*width];
	// initialise un tableau pseudo-aléatoire
	for(int i = 0; i < 8*width; i++) {
		A[i] = reduce(5*i ^ 17*i ^ 42);
	}

	// calcule les [width] FFTs parallèles en O(width * N^2).
	// B[i] = sum(A[j] * (omega^i)^j, i=0..127)
	for(int w=0; w<width; w++) {
		i16 omega_i = 1; // contient omega^i
		for(int i = 0; i < 8; i++) {
			i16 omega_ij = 1; // contient omega^(ij)
			B[w + i * width] = 0;
			for(int j = 0; j < 8; j++) {
				B[w + i * width] = reduce(B[w + i * width] + A[w + j*width] * omega_ij);
				omega_ij = reduce(omega_ij * omega_i);
			}
			omega_i = reduce(omega_i * omega);
		}
	}

        
	// check
	dif_fft8(A);
        
	// Output is in revbin order 
	for(int i=0; i<8; i++) {
	  for(int j=0; j<width; j++){
	    if (A[j + revbin8[i]*width] != B[j + i*width]) {
  printf("i : %d, j :%d a :%d, b : %d\n", i, j, A[j + revbin8[i]*width], B[j + i*width]); 
	      return 0;
	    }
	  }
	}
	return 1;

}


int main() {
  printf("parallel_reduce : %d\n", test_parallelreduce());
  printf("dif_butterfly :%d\n", test_dif_butterfly());
  printf("dif_fft8 : %d\n", test_fft8(8,4));
  return 0; 
}
