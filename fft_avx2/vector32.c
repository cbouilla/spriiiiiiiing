#include <stdlib.h>
#include <stdio.h>

//#include "compat.h"
#include "vector.h"

//#define PRINT_SOME 0

// #define CV(x) { .u16 = {x, x, x, x, x, x, x, x}}
 const v32 cst128 = v32_cst(128);
 const v32 cst255 = v32_cst(255);
 const v32 cst257 = v32_cst(257);
//static const union cv8  V0 = CV(0);


/*
 * Reduce modulo 257; result is in [-127; 383]
 * REDUCE(x) := (x&255) - (x>>8)
 */
#define REDUCE(x)                               \
  ((x & cst255) - v32_shift_r(x, 8))

/*
 * Reduce from [-127; 383] to [-128; 128]
 * EXTRA_REDUCE_S(x) := x<=128 ? x : x-257
 */
#define EXTRA_REDUCE(x)                       \
  (x - (cst257 & v32_cmp_gt(x, cst128)))

/*
 * Reduce modulo 257; result is in [-128; 128]
 */
#define REDUCE_FULL(x)                        \
  EXTRA_REDUCE(REDUCE(x))


// u, v, #bits to shift
#define DIF_BUTTERFLY(_u,_v, n)         \
  do {                             \
    const v32 u = _u;              \
    const v32 v = _v;              \
    _u = u + v;                    \
    if (n)                         \
      _v = v32_shift_l(u - v, n);  \
    else                           \
      _v = u - v;                  \
  } while(0)

#define DIT_BUTTERFLY(_u, _v, n)   \
  do {                             \
    const v32 u = _u;              \
    const v32 v = (n) ? v32_shift_l(_v, n) : _v; \
    _u = u + v;                    \
    _v = u - v;                    \
  } while(0)

#define DIT_TROUBLESOME_BUTTERFLY(_x, n)            \
  do {                                              \
    const v16 u = _mm256_extracti128_si256(_x, 0);  \
    const v16 v = _mm256_extracti128_si256(_x, 1);  \
    const v16 w = n ? v16_shift_l(v, n) : v;        \
    const v16 foo = u + w;                          \
    const v16 bar = u - w;                          \
    const v32 fooooo = _mm256_castsi128_si256(foo); \
    _x = _mm256_inserti128_si256 (fooooo, bar, 1);  \
  } while(0)

#define INTERLEAVE(x,y)                          \
  do {                                           \
    const v32 t1= x;                             \
    const v32 t2= y;                             \
    x = _mm256_unpacklo_epi16(t1, t2);           \
    y = _mm256_unpackhi_epi16(t1, t2);           \
  } while(0)

#define BIG_INTERLEAVE(_x,_y)                       \
  do {                                              \
    v32 t1= _x;                                      \
    v32 t2= _y;                                      \
    _x = _mm256_permute2x128_si256(t1, t2, 0x20);    \
    _y = _mm256_permute2x128_si256(t1, t2, 0x31);    \
  } while(0)


/* computes 16 parallel size-8 NTT with omega=4 */
// output in revbin order. GCC deals decently with this.
void dif_fft8(void *a) {
  v32* const A = a;
  register v32 X0, X1, X2, X3, X4, X5, X6, X7;

  X0 = A[0];
  X1 = A[1];
  X2 = A[2];
  X3 = A[3];
  X4 = A[4];
  X5 = A[5];
  X6 = A[6];
  X7 = A[7];

  // TODO : check bounds everywhere

  // outer layer
  DIF_BUTTERFLY(X0, X4, 0);
  DIF_BUTTERFLY(X1, X5, 2);
  DIF_BUTTERFLY(X2, X6, 4);
  DIF_BUTTERFLY(X3, X7, 6);
  
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);
  
  // middle layer
  DIF_BUTTERFLY(X0, X2, 0);
  DIF_BUTTERFLY(X4, X6, 0);
  DIF_BUTTERFLY(X1, X3, 4);
  DIF_BUTTERFLY(X5, X7, 4);
  
  // external layer
  DIF_BUTTERFLY(X0, X1, 0);
  DIF_BUTTERFLY(X2, X3, 0);
  DIF_BUTTERFLY(X4, X5, 0);
  DIF_BUTTERFLY(X6, X7, 0);
  
  /* We don't need to reduce X(0) */
  A[0] = REDUCE_FULL(X0); // FIXME later
  A[1] = REDUCE_FULL(X1);
  A[2] = REDUCE_FULL(X2);
  A[3] = REDUCE_FULL(X3);
  A[4] = REDUCE_FULL(X4);
  A[5] = REDUCE_FULL(X5);
  A[6] = REDUCE_FULL(X6);
  A[7] = REDUCE_FULL(X7);
}



/* computes 8 parallel size-16 NTT with omega=2 */
// input in revbin order.
void dit_fft16(void *a) {
  v32* A = a;
  register v32 X0, X1, X2, X3, X4, X5, X6, X7;
  
  X0 = A[0]; // 0-1
  X1 = A[1]; // 2-3 
  X2 = A[2]; // 4-5
  X3 = A[3]; // 6-7
  X4 = A[4]; // 8-9
  X5 = A[5]; // a-b
  X6 = A[6]; // c-d
  X7 = A[7]; // e-f

  BIG_INTERLEAVE(X0, X4);
  BIG_INTERLEAVE(X1, X5);
  BIG_INTERLEAVE(X2, X6);
  BIG_INTERLEAVE(X3, X7);

  // maintenant
  // X0 = 0, 8
  // X1 = 2, a
  // X2 = 4, c
  // X3 = 6, e
  // X4 = 1, 9
  // X5 = 3, b
  // X6 = 5, d
  // X7 = 7, f

  // TODO : vérifier les bornes partout, là-dedans...

 /* first phase
  DIT-butterfly 0 <---> 1  [omega=1]       X0 <--> X4 
  DIT-butterfly 8 <---> 9  [omega=1]       
  DIT-butterfly 2 <---> 3  [omega=1]       X1 <--> X5
  DIT-butterfly a <---> b  [omega=1]     
  DIT-butterfly 4 <---> 5  [omega=1]       X2 <--> X6
  DIT-butterfly c <---> d  [omega=1]
  DIT-butterfly 6 <---> 7  [omega=1]       X3 <---> X7
  DIT-butterfly d <---> f  [omega=1]
*/
  DIT_BUTTERFLY(X0, X4, 0);
  DIT_BUTTERFLY(X1, X5, 0);
  DIT_BUTTERFLY(X2, X6, 0);
  DIT_BUTTERFLY(X3, X7, 0);
 

/* second phase
  DIT-butterfly 0  <---> 2  [omega=1]     X0 <---> X1
  DIT-butterfly 8  <---> 10 [omega=1]  
  DIT-butterfly 4  <---> 6  [omega=1]     X2 <---> X3
  DIT-butterfly 12 <---> 14 [omega=1]
  DIT-butterfly 1  <---> 3  [omega=16]    X4 <---> X5
  DIT-butterfly 9  <---> 11 [omega=16] 
  DIT-butterfly 5  <---> 7  [omega=16]    X6 <---> X7
  DIT-butterfly 13 <---> 15 [omega=16]
  */
  DIT_BUTTERFLY(X0, X1, 0);
  DIT_BUTTERFLY(X2, X3, 0);
  DIT_BUTTERFLY(X4, X5, 4);
  DIT_BUTTERFLY(X6, X7, 4);

 // TODO : vérifier à quel point ceci est nécessaire
  X4 = REDUCE(X4);
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);

/* third phase
  #DIT-butterfly 0  <---> 4   [omega=1]     X0 <---> X2
  #DIT-butterfly 8  <---> 12  [omega=1]  
  #DIT-butterfly 1  <---> 5   [omega=4]     X4 <---> X6
  #DIT-butterfly 9  <---> 13  [omega=4]  
  #DIT-butterfly 2  <---> 6   [omega=16]    X1 <---> X3
  #DIT-butterfly 10 <---> 14  [omega=16]
  #DIT-butterfly 3  <---> 7   [omega=64]    X5 <---> X7
  #DIT-butterfly 11 <---> 15  [omega=64]
*/
  DIT_BUTTERFLY(X0, X2, 0);
  DIT_BUTTERFLY(X4, X6, 2);
  DIT_BUTTERFLY(X1, X3, 4);
  DIT_BUTTERFLY(X5, X7, 6);

 // TODO : vérifier à quel point ceci est nécessaire
  X4 = REDUCE(X4);
  X5 = REDUCE(X5);
  X1 = REDUCE(X1);
  X3 = REDUCE(X3);
  X5 = REDUCE(X5);
  X7 = REDUCE(X7);

  /* fourth phase
  DIT-butterfly 0 <---> 8   [omega=1]     
  DIT-butterfly 1 <---> 9   [omega=2]     
  DIT-butterfly 2 <---> 10  [omega=4]    
  DIT-butterfly 3 <---> 11  [omega=8]    
  DIT-butterfly 4 <---> 12  [omega=16]   
  DIT-butterfly 5 <---> 13  [omega=32]   
  DIT-butterfly 6 <---> 14  [omega=64]   
  DIT-butterfly 7 <---> 15  [omega=128]  
*/
  DIT_TROUBLESOME_BUTTERFLY(X0, 0);
  DIT_TROUBLESOME_BUTTERFLY(X4, 1);
  DIT_TROUBLESOME_BUTTERFLY(X1, 2);
  DIT_TROUBLESOME_BUTTERFLY(X5, 3);
  DIT_TROUBLESOME_BUTTERFLY(X2, 4);
  DIT_TROUBLESOME_BUTTERFLY(X6, 5);
  DIT_TROUBLESOME_BUTTERFLY(X3, 6);
  DIT_TROUBLESOME_BUTTERFLY(X7, 7);

  // remet dans le bon ordre
  BIG_INTERLEAVE(X0, X4);
  BIG_INTERLEAVE(X1, X5);
  BIG_INTERLEAVE(X2, X6);
  BIG_INTERLEAVE(X3, X7);

  /* We don't need to reduce X(0) */
  A[0] = REDUCE_FULL(X0);
  A[1] = REDUCE_FULL(X1);
  A[2] = REDUCE_FULL(X2);
  A[3] = REDUCE_FULL(X3);
  A[4] = REDUCE_FULL(X4);
  A[5] = REDUCE_FULL(X5);
  A[6] = REDUCE_FULL(X6);
  A[7] = REDUCE_FULL(X7);
}

static const v32 FFT128_Twiddle[8] =  {
  { 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1}, 
  { 1,  -60,    2, -120,    4,   17,    8,   34,   16,   68,   32, -121,   64,   15,  128,   30}, 
  { 1,  -35,  -60,   44,    2,  -70, -120,   88,    4,  117,   17,  -81,    8,  -23,   34,   95}, 
  { 1,   44, -120,  117,    8,   95,   68,  -92,   64,  -11,   30,   35,   -2,  -88,  -17,   23}, 
  { 1,   42,  -35,   72,  -60,   50,   44,   49,    2,   84,  -70, -113, -120,  100,   88,   98}, 
  { 1,   50,  -70,   98,   17,   79,   95,  124,   32,   58,   73,   52,   30,  -42,  -44,  113}, 
  { 1,   72,   44,   84, -120,   98,  117,  -57,    8,   62,   95,  -99,   68,   13,  -92,   58}, 
  { 1,   49,   88,  -57,   34,  124,  -92,  118,  128,  104,  -44, -100,  -17,  -62,   46,  -59}};

// the whole she-bang. FFT-128 with omega=42 !
// input in [-128,383] (i.e. "REDUCEd" but not "EXTRA_REDUCEd")
void fft128(void *a) {
  v32* const A = a;
  register v32 X0, X1, X2, X3, X4, X5, X6, X7;

  X0 = A[0];
  X1 = A[1];
  X2 = A[2];
  X3 = A[3];
  X4 = A[4];
  X5 = A[5];
  X6 = A[6];
  X7 = A[7];

  // STEP 1 : 16x parallel DIF FFT-8 with omega=4
  DIF_BUTTERFLY(X0, X4, 0);
  DIF_BUTTERFLY(X1, X5, 2);
  DIF_BUTTERFLY(X2, X6, 4);
  DIF_BUTTERFLY(X3, X7, 6);
  X5 = REDUCE(X5);
  X7 = REDUCE(X7);

  DIF_BUTTERFLY(X0, X2, 0);
  DIF_BUTTERFLY(X4, X6, 0);
  DIF_BUTTERFLY(X1, X3, 4);
  DIF_BUTTERFLY(X5, X7, 4);
  
  DIF_BUTTERFLY(X0, X1, 0);
  DIF_BUTTERFLY(X2, X3, 0);
  DIF_BUTTERFLY(X4, X5, 0);
  DIF_BUTTERFLY(X6, X7, 0);
  
  // X0 = REDUCE_FULL(X0); 
  X1 = REDUCE_FULL(X1);
  X2 = REDUCE_FULL(X2);
  X3 = REDUCE_FULL(X3);
  X4 = REDUCE_FULL(X4);
  X5 = REDUCE_FULL(X5);
  X6 = REDUCE_FULL(X6);
  X7 = REDUCE_FULL(X7);


  // STEP 2 : multiply by twiddle factors
  X1 *= FFT128_Twiddle[1];
  X2 *= FFT128_Twiddle[2];
  X3 *= FFT128_Twiddle[3];
  X4 *= FFT128_Twiddle[4];
  X5 *= FFT128_Twiddle[5];
  X6 *= FFT128_Twiddle[6];
  X7 *= FFT128_Twiddle[7];

  // qu'est-ce qui est strictement nécessaire là-dedans ?
  X0 = REDUCE(X0);
  X1 = REDUCE(X1);
  X2 = REDUCE(X2);
  X3 = REDUCE(X3);
  X4 = REDUCE(X4);
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);

  // STEP 3 : (nearly complete) transpose
  INTERLEAVE(X0, X1);
  INTERLEAVE(X2, X3);
  INTERLEAVE(X4, X5);
  INTERLEAVE(X6, X7);

  INTERLEAVE(X0, X2);
  INTERLEAVE(X1, X3);
  INTERLEAVE(X4, X6);
  INTERLEAVE(X5, X7);

  INTERLEAVE(X0, X4);
  INTERLEAVE(X1, X5);
  INTERLEAVE(X2, X6);
  INTERLEAVE(X3, X7);

  // STEP 4 : 8x parallel DIT FFT-16 with omega=2
  BIG_INTERLEAVE(X0, X4);
  BIG_INTERLEAVE(X1, X5);
  BIG_INTERLEAVE(X2, X6);
  BIG_INTERLEAVE(X3, X7);

  DIT_BUTTERFLY(X0, X4, 0);
  DIT_BUTTERFLY(X1, X5, 0);
  DIT_BUTTERFLY(X2, X6, 0);
  DIT_BUTTERFLY(X3, X7, 0);
 
  DIT_BUTTERFLY(X0, X1, 0);
  DIT_BUTTERFLY(X2, X3, 0);
  DIT_BUTTERFLY(X4, X5, 4);
  DIT_BUTTERFLY(X6, X7, 4);
  X4 = REDUCE(X4);
  X5 = REDUCE(X5);
  X6 = REDUCE(X6);
  X7 = REDUCE(X7);

  DIT_BUTTERFLY(X0, X2, 0);
  DIT_BUTTERFLY(X4, X6, 2);
  DIT_BUTTERFLY(X1, X3, 4);
  DIT_BUTTERFLY(X5, X7, 6);
  X2 = REDUCE(X2);
  X5 = REDUCE(X5);
  X1 = REDUCE(X1);
  X3 = REDUCE(X3);
  X5 = REDUCE(X5);
  X7 = REDUCE_FULL(X7);

// TODO : réfléchir à combiner BUTTERFLY et REDUCE

  DIT_TROUBLESOME_BUTTERFLY(X0, 0);
  DIT_TROUBLESOME_BUTTERFLY(X4, 1);
  DIT_TROUBLESOME_BUTTERFLY(X1, 2);
  DIT_TROUBLESOME_BUTTERFLY(X5, 3);
  DIT_TROUBLESOME_BUTTERFLY(X2, 4);
  DIT_TROUBLESOME_BUTTERFLY(X6, 5);
  DIT_TROUBLESOME_BUTTERFLY(X3, 6);
  DIT_TROUBLESOME_BUTTERFLY(X7, 7);

  // complete transpose
  BIG_INTERLEAVE(X0, X4);
  BIG_INTERLEAVE(X1, X5);
  BIG_INTERLEAVE(X2, X6);
  BIG_INTERLEAVE(X3, X7);

  /* We don't need to reduce X(0) */
  A[0] = REDUCE_FULL(X0);
  A[1] = REDUCE_FULL(X1);
  A[2] = REDUCE_FULL(X2);
  A[3] = REDUCE_FULL(X3);
  A[4] = REDUCE_FULL(X4);
  A[5] = REDUCE_FULL(X5);
  A[6] = REDUCE_FULL(X6);
  A[7] = REDUCE_FULL(X7);
} 

const v32 ZERO_VECT = v32_cst(0);
const v32 REJECTION_VECT = v32_cst(-1);
const v32 rm0 = {0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00};
const v32 rm1 = {0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0, 0x00, 0xf0};
const __v32qi p0 = {0x00, 0x04, 0x08, 0x0c, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
                    0x00, 0x04, 0x08, 0x0c, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};
const __v8si  p1 = {0, 4, 1, 2, 3, 5, 6, 7};
const __v32qi p2 = {0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00, 0x02, 0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e,
                    0x00, 0x02, 0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};


void dump(const char *name, __m256i a) {
  printf("%s = ", name);
  for(int i = 0; i < 32; i++) {
    printf("%02x-", (unsigned char) ((__v32qi) a)[i]);
  }
  printf("\n");  
}


void dump16(const char *name, __m128i a) {
  printf("%s = ", name);
  for(int i = 0; i < 16; i++) {
    printf("%02x-", (unsigned char) ((__v16qi) a)[i]);
  }
  printf("\n");  
}

// rejection sampling
int keep(const v32 a) {
  return v32_movemask(v32_cmp_eq(a, REJECTION_VECT));
}

// extract just the sign
uint16_t msb(const v32 a) {
  v32 mask = v32_cmp_gt(a, ZERO_VECT);
  v32 stacked = _mm256_shuffle_epi8(mask, p2);
  return v32_movemask(stacked) >> 8;
}

uint64_t rounding(const v32 a) {
  v32 b = _mm256_srli_epi32(a & rm0, 0);
  v32 c = _mm256_srli_epi32(a & rm1, 20);
  v32 d = b ^ c; // in d, each dword contains a meaningful top byte
  v32 e = _mm256_shuffle_epi8(d, p0); // inside each lane, everything is OK
  v32 f = _mm256_permutevar8x32_epi32(e, p1); // move everything to the first lane
  return _mm256_extract_epi64(f, 0); // get it
}

/*
 * multiply 128 bits (actually only first 127) with the generator matrix of bch code
 * in order to produce only 64 bits. the code distance is d=21.
 * we then use the 128th bit as a parity extension bit to the code and get distance d=22
 * input: array of 2 packed 64bit integers
 * output: packed 64 bit integer
 *
 * we have 4 possible generator polynomials:
 * [0, 2, 7, 8, 10, 12, 14, 15, 16, 23, 25, 27, 28, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 44, 45, 48, 58, 61, 63]
 * [0, 2, 5, 15, 18, 19, 21, 22, 23, 24, 25, 26, 30, 31, 32, 33, 35, 36, 38, 40, 47, 48, 49, 51, 53, 55, 56, 61, 63]
 * [0, 1, 2, 3, 5, 8, 13, 17, 19, 21, 23, 27, 28, 32, 34, 35, 36, 39, 41, 43, 44, 50, 52, 54, 59, 60, 61, 62, 63]
 * [0, 1, 2, 3, 4, 9, 11, 13, 19, 20, 22, 24, 27, 28, 29, 31, 35, 36, 40, 42, 44, 46, 50, 55, 58, 60, 61, 62, 63]
 *
 * no specific polynomial is preferable. so we choose the first one
 */
 const __v2di generatorBCHMatrix = {0xa40137e3da81d585, 0xa40137e3da81d585};
uint64_t BCH128to64_clmul (const __v2di in) {
  uint64_t parity = _mm_extract_epi64(in, 0);
  __m128i t0 = _mm_clmulepi64_si128(in, generatorBCHMatrix, 0x00);
  __m128i t1 = _mm_clmulepi64_si128(in, generatorBCHMatrix, 0x11);
  __m128i p = t0 ^ _mm_shuffle_epi32(t1,  0x4e);
  uint64_t low = _mm_extract_epi64(p, 0);
  return low ^ -(parity & 1);
}

uint64_t BCH128to64 (const __v2di in) {
  register uint64_t b1 = in[0];
  register uint64_t res = b1;
  res ^= (b1 << 2);
  res ^= (b1 << 7);
  res ^= (b1 << 8);
  res ^= (b1 << 10);
  res ^= (b1 << 12);
  res ^= (b1 << 14);
  res ^= (b1 << 15);
  res ^= (b1 << 16);
  res ^= (b1 << 23);
  res ^= (b1 << 25);
  res ^= (b1 << 27);
  res ^= (b1 << 28);
  res ^= (b1 << 30);
  res ^= (b1 << 31);
  res ^= (b1 << 32);
  res ^= (b1 << 33);
  res ^= (b1 << 37);
  res ^= (b1 << 38);
  res ^= (b1 << 39);
  res ^= (b1 << 40);
  res ^= (b1 << 41);
  res ^= (b1 << 42);
  res ^= (b1 << 44);
  res ^= (b1 << 45);
  res ^= (b1 << 48);
  res ^= (b1 << 58);
  res ^= (b1 << 61);
  res ^= (b1 << 63);

  register uint64_t b2 = in[1];
  res ^= (b2 >> 62);
  res ^= (b2 >> 57);
  res ^= (b2 >> 56);
  res ^= (b2 >> 54);
  res ^= (b2 >> 52);
  res ^= (b2 >> 50);
  res ^= (b2 >> 49);
  res ^= (b2 >> 48);
  res ^= (b2 >> 41);
  res ^= (b2 >> 39);
  res ^= (b2 >> 37);
  res ^= (b2 >> 36);
  res ^= (b2 >> 34);
  res ^= (b2 >> 33);
  res ^= (b2 >> 32);
  res ^= (b2 >> 31);
  res ^= (b2 >> 27);
  res ^= (b2 >> 26);
  res ^= (b2 >> 25);
  res ^= (b2 >> 24);
  res ^= (b2 >> 23);
  res ^= (b2 >> 22);
  res ^= (b2 >> 20);
  res ^= (b2 >> 19);
  res ^= (b2 >> 16);
  res ^= (b2 >> 6);
  res ^= (b2 >> 3);
  res ^= (b2 >> 1);

  return res ^ ((uint64_t)(-(b2&1)));
}