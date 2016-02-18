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


// Test Butterfly :
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

int test_dit_butterfly() {
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

 
  DIT_BUTTERFLY(a, b, 1);

 
  //simulate the butterfly

  for(int i = 0; i < 8; i++) {
    tmp = c[i];
    c[i] = c[i] + (d[i] << 2);
    d[i] = tmp - (d[i] << 2);
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

int test_fft(int N, int width, i16 omega) {
  i16 A[N*width] __attribute ((aligned (16)));
  i16 B[N*width];
	// initialise un tableau pseudo-aléatoire
	for(int i = 0; i < N*width; i++) {
		A[i] = reduce(5*i ^ 17*i ^ 42);
	}

	// calcule les [width] FFTs parallèles en O(width * N^2).
	// B[i] = sum(A[j] * (omega^i)^j, i=0..N-1)
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
	if(N == 8) {
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

	} else if(N == 64) {
	  fft64(A);

	  /* for(int i = 0; i < 64; i++){ */
	  /*   printf("A[%d] = %d B = %d\n", i, A[i], B[i]); */
	  /* } */

	  for(int i = 0; i < 64; i++){
	    if(A[i] != B[i]) {
	      printf("A[%d] = %d vs B[%d] = %d\n", i, A[i], i, B[i]);

	      return 0;
	    }
	  }
	  return 1;
	} else if (N == 128){
	  fft128(A);

	  for(int i = 0; i < 128; i++){
	    if(A[i] != B[i]) {
	      printf("A[%d] = %d vs B[%d] = %d\n", i, A[i], i, B[i]);
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

// Teste le mode compteur pour x = 0x2 et Gray = 3+1 = 4.
int test_UpdateGray() {
  uint32_t x = 0;

  x = UpdateGray(0x2, 4);

  if((x ^ 0x6) != 0){
    printf("x : %x\n", x);
    return 0;
  }

  return 1;

}


// Teste la fonction reject
int test_reject(v16 a){
  int r, check;

  r = reject(a);
  check = 0;

  for(int i = 0; i < 8; i++){
    if(a[i] == -1){
      check = 1;
    }
  }

  if(check == 0 && r != 0){
    printf("r : %x\n", r);
    return 0;
  }

  if(check != 0 && r == 0){
    printf("check = %d\n", check);
    return 0;
  }

  return 1;
}

int test_rounding(){
  v16 a;
  uint32_t r = 0; 

  //initialise un tableau pseudo-aléatoire
  for(int i = 0; i < 8; i++){
    a[i] = rand() & 0x00ff;
  }
  for(int i = 0; i < 8; i++) {
    r <<=4;
    r ^= (a[i]>>4);
  }
 
  uint64_t r2 = rounding4(a);
  return r == r2;
}

int main() {
  v16 a = CV(1);
  a[2] = -1;

  printf("parallel_reduce : %d\n", test_parallelreduce());
  printf("dif_butterfly : %d\n", test_dif_butterfly());
  printf("dit_butterfly : %d\n", test_dit_butterfly());
  printf("dif_fft8 : %d\n", test_fft(8,8,4));
  printf("fft64 : %d\n", test_fft(64, 1, 46)); 
  printf("fft128 : %d\n", test_fft(128, 1, -118));
  printf("UpdateGray : %d\n", test_UpdateGray());
  printf("test_reject : %d\n", test_reject(a));
  printf("test_rounding : %d\n", test_rounding());
  return 0; 
}
