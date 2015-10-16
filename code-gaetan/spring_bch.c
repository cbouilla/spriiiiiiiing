#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include "spring.h"
#include "vector.c"
#include "rand_generated_arrays.c"

typedef v64 PackedBits;
typedef v64 R2poly;


// some constant definitions

// size of Z_q
#define FIELD_SIZE 257

// a generator for the Z_q field
#define GENERATOR 3
// discrete log of inverse(N) in the field (the base of the exponent is GENERATOR)
#define LOG_OF_INV_N 176
// the chosen omega for the DFT.
// used to get the principal evaluations on omega**t where t is odd
// of order 2*N
//#define OMEGA 82
#define OMEGA 41

#if N != 128
#error Only N=128 is supported
#endif

#if ((K!=64) && (K!=128))
#error Only K in {64,128} is supported
#endif


// a table of the powers of GENERATOR.	i --> GENERATOR ** I	 ( mod FIELD_SIZE )
static const int16_t generatorPowers[257] = {
	1, 3, 9, 27, 81, -14, -42, -126, -121, -106, -61, 74, -35, -105,
	-58, 83, -8, -24, -72, 41, 123, 112, 79, -20, -60, 77, -26, -78,
	23, 69, -50, 107, 64, -65, 62, -71, 44, -125, -118, -97, -34,
	-102, -49, 110, 73, -38, -114, -85, 2, 6, 18, 54, -95, -28, -84,
	5, 15, 45, -122, -109, -70, 47, -116, -91, -16, -48, 113, 82,
	-11, -33, -99, -40, -120, -103, -52, 101, 46, -119, -100, -43,
	128, 127, 124, 115, 88, 7, 21, 63, -68, 53, -98, -37, -111, -76,
	29, 87, 4, 12, 36, 108, 67, -56, 89, 10, 30, 90, 13, 39, 117, 94,
	25, 75, -32, -96, -31, -93, -22, -66, 59, -80, 17, 51, -104, -55,
	92, 19, 57, -86, -1, -3, -9, -27, -81, 14, 42, 126, 121, 106, 61,
	-74, 35, 105, 58, -83, 8, 24, 72, -41, -123, -112, -79, 20, 60,
	-77, 26, 78, -23, -69, 50, -107, -64, 65, -62, 71, -44, 125, 118,
	97, 34, 102, 49, -110, -73, 38, 114, 85, -2, -6, -18, -54, 95,
	28, 84, -5, -15, -45, 122, 109, 70, -47, 116, 91, 16, 48, -113,
	-82, 11, 33, 99, 40, 120, 103, 52, -101, -46, 119, 100, 43, -128,
	-127, -124, -115, -88, -7, -21, -63, 68, -53, 98, 37, 111, 76,
	-29, -87, -4, -12, -36, -108, -67, 56, -89, -10, -30, -90, -13,
	-39, -117, -94, -25, -75, 32, 96, 31, 93, 22, 66, -59, 80, -17,
	-51, 104, 55, -92, -19, -57, 86, 1};
static const v8 generatorPowersT1 = 
{0x01, 0x03, 0x09, 0x1b, 0x51, 0xf2, 0xd6, 0x82,
 0x87, 0x96, 0xc3, 0x4a, 0xdd, 0x97, 0xc6, 0x53};
static const v8 generatorPowersT2 =
{0xff, 0xf6, 0x3e, 0x0, 0xee, 0x7e, 0x02, 0xde,
 0xfd, 0x6, 0xbe, 0xfc, 0xe, 0x7f, 0xfa, 0x1e};


// a table of the powers of omega.	i --> omega**(2N-i)		( mod FIELD_SIZE )

static const v16 omegaPowers[] = {
	{1, -94, 98, 40, 95, 65, 58, -55},
	{30, 7, 113, -85, 23, -106, -59, -108},
	{-128, -47, 49, 20, -81, -96, 29, 101},
	{15, -125, -72, 86, -117, -53, 99, -54},
	{-64, 105, -104, 10, 88, -48, -114, -78},
	{-121, 66, -36, 43, 70, 102, -79, -27},
	{-32, -76, -52, 5, 44, -24, -57, -39},
	{68, 33, -18, -107, 35, 51, 89, 115},
	{-16, -38, -26, -126, 22, -12, 100, 109},
	{34, -112, -9, 75, -111, -103, -84, -71},
	{-8, -19, -13, -63, 11, -6, 50, -74},
	{17, -56, 124, -91, 73, 77, -42, 93},
	{-4, 119, 122, 97, -123, -3, 25, -37},
	{-120, -28, 62, 83, -92, -90, -21, -82},
	{-2, -69, 61, -80, 67, 127, -116, 110},
	{-60, -14, 31, -87, -46, -45, 118, -41}
};

// tables containing the exponents which represent the polynomials (our function key)
uint8_t A[K+1][N] ALIGN;
v8 AA[K/8][256][N/16] ALIGN;
uint8_t B[K+1][N/2] ALIGN;
int16_t D[K+1][2][N] ALIGN; // D[k][1] is inv

R2poly multiplyR2PolynomialsInRing (R2poly a, R2poly b)
{
#if USE_CLMUL

	R2poly aa = v64_xor((v64)v32_shufxor((v32)a,2),a);
	// aa = {a[0]^a[1], a[0]^a[1]};
	R2poly bb = v64_xor((v64)v32_shufxor((v32)b,2),b);
	// bb = {b[0]^b[1], b[0]^b[1]};

	R2poly t0 = v64_clmul(a,b,0x00);
	R2poly t1 = v64_clmul(a,b,0x11);
	R2poly t2 = v64_clmul(aa,bb,0x00);

	R2poly t = v64_xor(t0,t1);
	t = v64_xor((v64)v32_shufxor((v32)t,2),t);
	// t = { t[0]^t[1], t[0]^t[1] };
	return v64_xor((v64)v32_shufxor((v32)t2,2),t);
	// return { t2[1]^t[0], t2[0]^t[1] };

#else

	union u64 x = {.v = a}, y = {.v = b}, xy = {.u = {0, 0}};

	for (int i=0; i<128; i++) {
	  if (y.u[0]&1)
	    xy.v = v64_xor(xy.v, x.v);

		union u64 t = x;
	  x.u[0] = t.u[0]<<1 | t.u[1]>>63;
	  x.u[1] = t.u[1]<<1 | t.u[0]>>63;

		y.u[0] = y.u[0]>>1 | y.u[1]<<63;
		y.u[1] = y.u[1]>>1;

	}

	return xy.v;

#endif
}

#if USE_CLMUL == 0
R2poly multiplyR2PolynomialsFromTable (R2poly a, R2poly mulTables[256]) {
	union cv8 t = {.v8 = (v8)a};
	R2poly mul;

	R2poly res = mulTables[t.u8[0]];

#define MUL(i) do {																\
		mul = mulTables[t.u8[16-i]];									\
		mul = v64_rotate_bytes(mul,i);								\
		res = v64_xor (res, mul);											\
	} while(0)

	MUL( 1); 	MUL( 2); 	MUL( 3); 	MUL( 4);
	MUL( 5); 	MUL( 6); 	MUL( 7); 	MUL( 8);
	MUL( 9); 	MUL(10); 	MUL(11); 	MUL(12);
	MUL(13); 	MUL(14); 	MUL(15);

	return res;
}
#endif

int degree (R2poly a)
{
	union u64 u = { .v= a};
  return u.u[1]? 127-__builtin_clzll(u.u[1]):
    u.u[0]? 63-__builtin_clzll(u.u[0]): -1;
}

void printPolynomial (R2poly a)
{
	union u64 u = { .v=a };
  printf("0x%016" PRIx64 "%016" PRIx64, u.u[1], u.u[0]);
}

int parity (R2poly a)
{
	union u64 u = { .v= a};
	return __builtin_parityll(u.u[0]^u.u[1]);
}

int isInvertible (R2poly a)
{
	return parity(a) == 1;
}


R2poly cyclicShiftLeft (R2poly a, int count)
{
	count %= N;
	union u64 x = {.v = a};
	
	if (count >= 64) {
		uint64_t tmp = x.u[1];
		x.u[1] = x.u[0];
		x.u[0] = tmp;
		count -= 64;
	}
	if (count != 0)	{
		uint64_t tmp = x.u[1];
		x.u[1] <<= count;
		x.u[1] |= (x.u[0] >> (64-count));
		x.u[0] <<= count;
		x.u[0] |= (tmp >> (64-count));
	}
	
	return x.v;
}

R2poly inverseInRing2 (R2poly poly)
{
	const union u64 one = {.u = {1, 0}};
	const union u64 zero = {.u = {0, 0}};
	assert (isInvertible(poly));

	int quotientShift = N - degree(poly);
	R2poly b = poly;
	R2poly a = cyclicShiftLeft(poly, quotientShift);
	R2poly y = one.v;
	R2poly lasty = cyclicShiftLeft (y, quotientShift);
	if (degree(a) < degree(b))
	{
		register R2poly tmp = a;
		a = b;
		b = tmp;
		tmp = y;
		y = lasty;
		lasty = tmp;
	}

	// Loop invariant: a*y  == b*lasty

	while (! v64_eq(b, zero.v)) {
		quotientShift = degree(a) - degree(b);
		a = v64_xor (a, cyclicShiftLeft(b, quotientShift));
		assert (degree(y) + quotientShift <= N);
		lasty = v64_xor (lasty, cyclicShiftLeft(y, quotientShift));

		if (degree(a) < degree(b))
		{
			// swap a, b
			// swap y, lasty
			register R2poly tmp = a;
			a = b;
			b = tmp;
			tmp = y;
			y = lasty;
			lasty = tmp;
		}
	}

	a = multiplyR2PolynomialsInRing(poly, lasty);
	assert (v64_eq(a, one.v));
	return lasty;
}


// initialize the new R_q polynomials representation table from table A
void initializeD ()
{
	for (int i=0; i<K+1; i++) {
		for (int j=0; j<N; j++) {
			D[i][0][j]		 = generatorPowers[A[i][j]];
			D[i][1][j] = generatorPowers[(256-A[i][j])%256];
		}
	}
}


/*
 * next function switches from a R_2 exponent representation of an invertible polynomial
 * into the radix base ( (1+x)**j for j=0,1,... ) representation.
 * This function is specifically only for the N=128 case.
 * input: expsArr[] - an array of the exponents of some invertible polynomial in R_2.
 *										the exponents are of the generators of the form 1+(1+x)**t where t is odd
 *										every element in this array is an exponent of one of this generators.
 *										length of the array is 64 (N/2)
 * output: packedRadixCode[] - an array of two uint64_t. I.e. 128 bits representing the
 *														 coefficients of the same polynomial in the radix base
 * the method is presented in algorithm 4.1 in the paper
 */
R2poly exponentsToRadixRepresentation (v8 subsetSumOfLogs[])
{
	uint64_t packedRadixCode[2] ALIGN;

	// radixBits <- e_1
	packedRadixCode[0] = 1;

#define TABLE1(z,a)              z^0,                     z^(1ULL<<a)
#define TABLE2(z,a,x)            TABLE1(z,a),             TABLE1(z^(1ULL<<x),a)
#define TABLE3(z,a,b,x)          TABLE2(z,a,b),           TABLE2(z^(1ULL<<x),a,b)
#define TABLE4(z,a,b,c,x)        TABLE3(z,a,b,c),         TABLE3(z^(1ULL<<x),a,b,c)
#define TABLE5(z,a,b,c,d,x)      TABLE4(z,a,b,c,d),       TABLE4(z^(1ULL<<x),a,b,c,d)
#define TABLE6(z,a,b,c,d,e,x)    TABLE5(z,a,b,c,d,e),     TABLE5(z^(1ULL<<x),a,b,c,d,e)
#define TABLE7(z,a,b,c,d,e,f,x)  TABLE6(z,a,b,c,d,e,f),   TABLE6(z^(1ULL<<x),a,b,c,d,e,f)
#define TABLE8(a,b,c,d,e,f,g,x)  TABLE7(0,a,b,c,d,e,f,g), TABLE7(  (1ULL<<x),a,b,c,d,e,f,g)

	// these code lines implement the algorithm lines:
	//		 for every valid (i,k):
	//				 if e_{i,k,lg(N)-i-1} == 1 then radixBits[N/2 + kN/2**(i+1)] <- 1

	uint8_t *expsArr = (uint8_t*) subsetSumOfLogs;

#if EXP_TABLES

	uint16_t idx = v8_extract_bits(subsetSumOfLogs[0], 6, 5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2);
	static const uint64_t T1[] = {TABLE8( 0, 32, 16, 48,  8, 24, 40, 56)};
	static const uint64_t T2[] = {TABLE8( 4, 12, 20, 28, 36, 44, 52, 60)};

	packedRadixCode[1]  = T1[idx&255];
	packedRadixCode[1] ^= T2[idx>>8];

	idx = v8_extract_bit(subsetSumOfLogs[1], 1);
	static const uint64_t T3[] = {TABLE8( 2,  6, 10, 14, 18, 22, 26, 30)};
	static const uint64_t T4[] = {TABLE8(34, 38, 42, 46, 50, 54, 58, 62)};

	packedRadixCode[1] ^= T3[idx&255];
	packedRadixCode[1] ^= T4[idx>>8];

	idx = v8_extract_bit(subsetSumOfLogs[2], 0);
	static const uint64_t T5[] = {TABLE8( 1,  3,  5,  7,  9, 11, 13, 15)};
	static const uint64_t T6[] = {TABLE8(17, 19, 21, 23, 25, 27, 29, 31)};

	packedRadixCode[1] ^= T5[idx&255];
	packedRadixCode[1] ^= T6[idx>>8];

	idx = v8_extract_bit(subsetSumOfLogs[3], 0);
	static const uint64_t T7[] = {TABLE8(33, 35, 37, 39, 41, 43, 45, 47)};
	static const uint64_t T8[] = {TABLE8(49, 51, 53, 55, 57, 59, 61, 63)};

	packedRadixCode[1] ^= T7[idx&255];
	packedRadixCode[1] ^= T8[idx>>8];

#else
	packedRadixCode[1] = (expsArr[0] & 64) >> 6;
	packedRadixCode[1] ^=	 ((uint64_t)(expsArr[1] & 32) << 27)
		^	 ((uint64_t)(expsArr[2] & 16) << 12)
		^	 ((uint64_t)(expsArr[3] & 16) << 44)
		^	 ((uint64_t)(expsArr[4] & 8) << 5)
		^	 ((uint64_t)(expsArr[5] & 8) << 21)
		^	 ((uint64_t)(expsArr[6] & 8) << 37)
		^	 ((uint64_t)(expsArr[7] & 8) << 53)
		^	 ((uint64_t)(expsArr[8] & 4) << 2)
		^	 ((uint64_t)(expsArr[9] & 4) << 10)
		^	 ((uint64_t)(expsArr[10] & 4) << 18)
		^	 ((uint64_t)(expsArr[11] & 4) << 26)
		^	 ((uint64_t)(expsArr[12] & 4) << 34)
		^	 ((uint64_t)(expsArr[13] & 4) << 42)
		^	 ((uint64_t)(expsArr[14] & 4) << 50)
		^	 ((uint64_t)(expsArr[15] & 4) << 58)
		^	 ((uint64_t)(expsArr[16] & 2) << 1)
		^	 ((uint64_t)(expsArr[17] & 2) << 5)
		^	 ((uint64_t)(expsArr[18] & 2) << 9)
		^	 ((uint64_t)(expsArr[19] & 2) << 13)
		^	 ((uint64_t)(expsArr[20] & 2) << 17)
		^	 ((uint64_t)(expsArr[21] & 2) << 21)
		^	 ((uint64_t)(expsArr[22] & 2) << 25)
		^	 ((uint64_t)(expsArr[23] & 2) << 29)
		^	 ((uint64_t)(expsArr[24] & 2) << 33)
		^	 ((uint64_t)(expsArr[25] & 2) << 37)
		^	 ((uint64_t)(expsArr[26] & 2) << 41)
		^	 ((uint64_t)(expsArr[27] & 2) << 45)
		^	 ((uint64_t)(expsArr[28] & 2) << 49)
		^	 ((uint64_t)(expsArr[29] & 2) << 53)
		^	 ((uint64_t)(expsArr[30] & 2) << 57)
		^	 ((uint64_t)(expsArr[31] & 2) << 61)
		^	 ((uint64_t)(expsArr[32] & 1) << 1)
		^	 ((uint64_t)(expsArr[33] & 1) << 3)
		^	 ((uint64_t)(expsArr[34] & 1) << 5)
		^	 ((uint64_t)(expsArr[35] & 1) << 7)
		^	 ((uint64_t)(expsArr[36] & 1) << 9)
		^	 ((uint64_t)(expsArr[37] & 1) << 11)
		^	 ((uint64_t)(expsArr[38] & 1) << 13)
		^	 ((uint64_t)(expsArr[39] & 1) << 15)
		^	 ((uint64_t)(expsArr[40] & 1) << 17)
		^	 ((uint64_t)(expsArr[41] & 1) << 19)
		^	 ((uint64_t)(expsArr[42] & 1) << 21)
		^	 ((uint64_t)(expsArr[43] & 1) << 23)
		^	 ((uint64_t)(expsArr[44] & 1) << 25)
		^	 ((uint64_t)(expsArr[45] & 1) << 27)
		^	 ((uint64_t)(expsArr[46] & 1) << 29)
		^	 ((uint64_t)(expsArr[47] & 1) << 31)
		^	 ((uint64_t)(expsArr[48] & 1) << 33)
		^	 ((uint64_t)(expsArr[49] & 1) << 35)
		^	 ((uint64_t)(expsArr[50] & 1) << 37)
		^	 ((uint64_t)(expsArr[51] & 1) << 39)
		^	 ((uint64_t)(expsArr[52] & 1) << 41)
		^	 ((uint64_t)(expsArr[53] & 1) << 43)
		^	 ((uint64_t)(expsArr[54] & 1) << 45)
		^	 ((uint64_t)(expsArr[55] & 1) << 47)
		^	 ((uint64_t)(expsArr[56] & 1) << 49)
		^	 ((uint64_t)(expsArr[57] & 1) << 51)
		^	 ((uint64_t)(expsArr[58] & 1) << 53)
		^	 ((uint64_t)(expsArr[59] & 1) << 55)
		^	 ((uint64_t)(expsArr[60] & 1) << 57)
		^	 ((uint64_t)(expsArr[61] & 1) << 59)
		^	 ((uint64_t)(expsArr[62] & 1) << 61)
		^	 ((uint64_t)(expsArr[63] & 1) << 63);
#endif

	// mask is used in order to avoid inefficient if-else clauses.
	// mask is always determined by 1 bit. if it is 0, then the bit-mask will be 000000000000....
	// if it is 1, then the bit-mask will be 11111111............
	// using bitwise & on some value along with the computed mask, determines if the
	// packedRadixCode[] values will be changed or not.
	register uint64_t mask;

	// these code lines implement the algorithm lines:
	//		 for every checked bit e_{i,k,l}:
	//				 if e_{i,k,l} == 1 then radixBits <- radixBits ^ (radixBits >> ((2**i+k)*2**l))
	mask = -(uint64_t)((expsArr[0] >> 5) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 32) | (packedRadixCode[0] >> 32)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 32) & mask);
	mask = -(uint64_t)((expsArr[0] >> 4) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 16) | (packedRadixCode[0] >> 48)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 16) & mask);
	mask = -(uint64_t)((expsArr[0] >> 3) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 8) | (packedRadixCode[0] >> 56)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 8) & mask);
	mask = -(uint64_t)((expsArr[0] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 4) | (packedRadixCode[0] >> 60)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 4) & mask);
	mask = -(uint64_t)((expsArr[0] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 2) | (packedRadixCode[0] >> 62)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 2) & mask);
	mask = -(uint64_t)(expsArr[0] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 1) | (packedRadixCode[0] >> 63)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 1) & mask);
	mask = -(uint64_t)((expsArr[1] >> 4) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 48) | (packedRadixCode[0] >> 16)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 48) & mask);
	mask = -(uint64_t)((expsArr[1] >> 3) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 24) | (packedRadixCode[0] >> 40)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 24) & mask);
	mask = -(uint64_t)((expsArr[1] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 12) | (packedRadixCode[0] >> 52)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 12) & mask);
	mask = -(uint64_t)((expsArr[1] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 6) | (packedRadixCode[0] >> 58)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 6) & mask);
	mask = -(uint64_t)(expsArr[1] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 3) | (packedRadixCode[0] >> 61)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 3) & mask);
	mask = -(uint64_t)((expsArr[2] >> 3) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 40) | (packedRadixCode[0] >> 24)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 40) & mask);
	mask = -(uint64_t)((expsArr[2] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 20) | (packedRadixCode[0] >> 44)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 20) & mask);
	mask = -(uint64_t)((expsArr[2] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 10) | (packedRadixCode[0] >> 54)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 10) & mask);
	mask = -(uint64_t)(expsArr[2] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 5) | (packedRadixCode[0] >> 59)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 5) & mask);
	mask = -(uint64_t)((expsArr[3] >> 3) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 56) | (packedRadixCode[0] >> 8)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 56) & mask);
	mask = -(uint64_t)((expsArr[3] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 28) | (packedRadixCode[0] >> 36)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 28) & mask);
	mask = -(uint64_t)((expsArr[3] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 14) | (packedRadixCode[0] >> 50)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 14) & mask);
	mask = -(uint64_t)(expsArr[3] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 7) | (packedRadixCode[0] >> 57)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 7) & mask);
	mask = -(uint64_t)((expsArr[4] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 36) | (packedRadixCode[0] >> 28)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 36) & mask);
	mask = -(uint64_t)((expsArr[4] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 18) | (packedRadixCode[0] >> 46)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 18) & mask);
	mask = -(uint64_t)(expsArr[4] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 9) | (packedRadixCode[0] >> 55)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 9) & mask);
	mask = -(uint64_t)((expsArr[5] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 44) | (packedRadixCode[0] >> 20)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 44) & mask);
	mask = -(uint64_t)((expsArr[5] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 22) | (packedRadixCode[0] >> 42)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 22) & mask);
	mask = -(uint64_t)(expsArr[5] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 11) | (packedRadixCode[0] >> 53)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 11) & mask);
	mask = -(uint64_t)((expsArr[6] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 52) | (packedRadixCode[0] >> 12)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 52) & mask);
	mask = -(uint64_t)((expsArr[6] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 26) | (packedRadixCode[0] >> 38)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 26) & mask);
	mask = -(uint64_t)(expsArr[6] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 13) | (packedRadixCode[0] >> 51)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 13) & mask);
	mask = -(uint64_t)((expsArr[7] >> 2) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 60) | (packedRadixCode[0] >> 4)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 60) & mask);
	mask = -(uint64_t)((expsArr[7] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 30) | (packedRadixCode[0] >> 34)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 30) & mask);
	mask = -(uint64_t)(expsArr[7] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 15) | (packedRadixCode[0] >> 49)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 15) & mask);
	mask = -(uint64_t)((expsArr[8] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 34) | (packedRadixCode[0] >> 30)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 34) & mask);
	mask = -(uint64_t)(expsArr[8] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 17) | (packedRadixCode[0] >> 47)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 17) & mask);
	mask = -(uint64_t)((expsArr[9] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 38) | (packedRadixCode[0] >> 26)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 38) & mask);
	mask = -(uint64_t)(expsArr[9] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 19) | (packedRadixCode[0] >> 45)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 19) & mask);
	mask = -(uint64_t)((expsArr[10] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 42) | (packedRadixCode[0] >> 22)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 42) & mask);
	mask = -(uint64_t)(expsArr[10] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 21) | (packedRadixCode[0] >> 43)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 21) & mask);
	mask = -(uint64_t)((expsArr[11] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 46) | (packedRadixCode[0] >> 18)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 46) & mask);
	mask = -(uint64_t)(expsArr[11] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 23) | (packedRadixCode[0] >> 41)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 23) & mask);
	mask = -(uint64_t)((expsArr[12] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 50) | (packedRadixCode[0] >> 14)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 50) & mask);
	mask = -(uint64_t)(expsArr[12] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 25) | (packedRadixCode[0] >> 39)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 25) & mask);
	mask = -(uint64_t)((expsArr[13] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 54) | (packedRadixCode[0] >> 10)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 54) & mask);
	mask = -(uint64_t)(expsArr[13] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 27) | (packedRadixCode[0] >> 37)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 27) & mask);
	mask = -(uint64_t)((expsArr[14] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 58) | (packedRadixCode[0] >> 6)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 58) & mask);
	mask = -(uint64_t)(expsArr[14] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 29) | (packedRadixCode[0] >> 35)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 29) & mask);
	mask = -(uint64_t)((expsArr[15] >> 1) & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 62) | (packedRadixCode[0] >> 2)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 62) & mask);
	mask = -(uint64_t)(expsArr[15] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 31) | (packedRadixCode[0] >> 33)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 31) & mask);
	mask = -(uint64_t)(expsArr[16] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 33) | (packedRadixCode[0] >> 31)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 33) & mask);
	mask = -(uint64_t)(expsArr[17] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 35) | (packedRadixCode[0] >> 29)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 35) & mask);
	mask = -(uint64_t)(expsArr[18] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 37) | (packedRadixCode[0] >> 27)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 37) & mask);
	mask = -(uint64_t)(expsArr[19] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 39) | (packedRadixCode[0] >> 25)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 39) & mask);
	mask = -(uint64_t)(expsArr[20] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 41) | (packedRadixCode[0] >> 23)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 41) & mask);
	mask = -(uint64_t)(expsArr[21] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 43) | (packedRadixCode[0] >> 21)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 43) & mask);
	mask = -(uint64_t)(expsArr[22] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 45) | (packedRadixCode[0] >> 19)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 45) & mask);
	mask = -(uint64_t)(expsArr[23] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 47) | (packedRadixCode[0] >> 17)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 47) & mask);
	mask = -(uint64_t)(expsArr[24] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 49) | (packedRadixCode[0] >> 15)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 49) & mask);
	mask = -(uint64_t)(expsArr[25] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 51) | (packedRadixCode[0] >> 13)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 51) & mask);
	mask = -(uint64_t)(expsArr[26] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 53) | (packedRadixCode[0] >> 11)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 53) & mask);
	mask = -(uint64_t)(expsArr[27] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 55) | (packedRadixCode[0] >> 9)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 55) & mask);
	mask = -(uint64_t)(expsArr[28] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 57) | (packedRadixCode[0] >> 7)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 57) & mask);
	mask = -(uint64_t)(expsArr[29] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 59) | (packedRadixCode[0] >> 5)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 59) & mask);
	mask = -(uint64_t)(expsArr[30] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 61) | (packedRadixCode[0] >> 3)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 61) & mask);
	mask = -(uint64_t)(expsArr[31] & 1);
	packedRadixCode[1] ^= (((packedRadixCode[1] << 63) | (packedRadixCode[0] >> 1)) & mask);
	packedRadixCode[0] ^= ((packedRadixCode[0] << 63) & mask);

	union u64 res = {.u = {packedRadixCode[0], packedRadixCode[1]}};
	return res.v;
}




/*
 * Input: subsetSumOfLogs[] - the summation of the generators exponents representing the R_2 polynomials
 *														(as given in table B). length of the array is N/2
 * output: R2poly - Packed R_2 polynomial bits - the coefficients of the matching
 *																					polynomial in power basis representation
 */
R2poly subsetSumOfExpsToR2Polynomial(v8 subsetSumOfLogs[])
{
	// change from exponent representation to radix basis representation
	union u64 representationPacked;
	representationPacked.v = exponentsToRadixRepresentation(subsetSumOfLogs);

	// change basis from radix basis ( 1+(1+x)**j for j=0,1,.....) representation
	// to power basis (x**j for j=0,1, ... ) representation
	// this is an iterative code for the recursive algorithm 4.2 in the paper
	const union u64 mask4 = { .u= {0xf0f0f0f0f0f0f0f0LL, 0xf0f0f0f0f0f0f0f0LL}};
	const union u64 mask2 = { .u= {0xccccccccccccccccLL, 0xccccccccccccccccLL}};
	const union u64 mask1 = { .u= {0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL}};

	representationPacked.u[0] ^= representationPacked.u[1];
	representationPacked.v ^=      v64_shift_r(     representationPacked.v, 32);
	representationPacked.v ^= (v64)v32_shift_r((v32)representationPacked.v, 16);
	representationPacked.v ^= (v64)v16u_shift_r((v16)representationPacked.v, 8);
	representationPacked.v ^= v64_shift_r(representationPacked.v & mask4.v, 4);
	representationPacked.v ^= v64_shift_r(representationPacked.v & mask2.v, 2);
	representationPacked.v ^= v64_shift_r(representationPacked.v & mask1.v, 1);

	return representationPacked.v;
}



// initialize the R_q polynomials exponent representation table from a (hard-coded) generated sequence of numbers
// this table is a part of the function key
void initializeA ()
{
	int i, j, r;

	r = 0;
	for (i = 0; i < K+1; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			A[i][j] = _randGeneratedNumbersForRq[r];
			++r;
		}
	}

	// A[K] represents the polynomial whose index is K.
	// This polynomial is ALWAYS present in the multiplication.
	// thus, this addition always adds exactly once LOG_OF_INV_N to the subset sum of logs.
	// this is the same as multiplying the result by 1/N which is an important factor in
	// the inverse DFT. We could just assume that this addition is already depicted in the
	// random generated numbers of the A table, and ignore the next addition.
	// Nevertheless it is done only once and does not influence performance. I decided to
	// add this, nevertheless, in order for the output to match our simulations.
	for (j = 0; j < N; ++j)
	{
		A[K][j] += LOG_OF_INV_N;
	}
}

/*
 * using (hard-coded) generated sequence of numbers
 * initialize table B which has the R_2 polynomials exponent representation
 * also initialize table C which has the R_2 polynomials themselves
 * in standard coefficient representation (and their inverses)
 * Although they represent the exact same polynomials, table C is used for faster
 * computation of the gray-code counter mode - but only if the chip supports
 * carry-less multiplication
 */
void initializeB_C ()
{
	int i, j, r;

	r = 0;
	for (i = 0; i < K+1; ++i)
	{
		for (j = 0; j < N/2; ++j)
		{
			B[i][j] = _randGeneratedNumbersForR2[r];
			++r;
		}


		R2poly poly = subsetSumOfExpsToR2Polynomial((v8*)B[i]);
		assert (isInvertible(poly));
	}

}

void initializeXX ()
{
	v8 (*vA)[N/16] = (void*)A;

	for (int kk=0; kk<K/8; kk++) {
		for (int i=0; i<N/16; i++) {
			AA[kk][0][i] ^= AA[kk][0][i];
			for (int j=0; j<8; j++) {
				AA[kk][1<<j][i] = vA[8*kk+j][i];
			}
			// Compute by linearity
			for (int k = 1; k<256; k++) {
				int msb = 1<<(__builtin_ffs(k)-1);
				assert(k&msb);
				AA[kk][k][i] = AA[kk][k-msb][i] + AA[kk][msb][i];
			}
		}
	}

	for (int i=0; i<N/16; i++) {
		for (int k=0; k<256; k++) {
			AA[0][k][i] += vA[K][i];
		}
	}
}


/*
 * multiplyRqPolynomials
 * This function returns the coefficients of the R_q multiplied polynomial
 * input of this function is the already calculate sum of logs of the evaluations of the
 * multiplier polynomials. It was chosen this way in order to enable fast computation of
 * the subset-sum, if we use the Gray-code counter mode.
 * Input - array of N uint8_t. the discrete logs of the result polynomial principal
 *				 evaluations on omega**t (when t is odd)
 * Output - the R_q polynomial coefficients
 */
void multiplyRqPolynomialsFromPoints(v16 output[], const v16 resRq[])
{
	for (int i=0; i<N/8; i++)
		output[i] = resRq[i];

	// change from the point-value evaluation space to the coefficient space
	//fftSIMDExplicitInstructions(output, input);
	fft128(output);

	// multiply by powers of OMEGA to get the actual coefficients of the polynomial multiplication
	// apply final modulo. The output is in canonical form, i.e. in the range 0..256
	// the added value CONGRUENCE_ZERO_I was chosen in such a way that the result will always
	// be positive but still it will not cause an overflow of the previous computations
	v16* vecOutputPtr = (v16*)output;

	for (int i = 0; i < N/8; ++i)
	{
		vecOutputPtr[i] = REDUCE(vecOutputPtr[i] * omegaPowers[i]);
	}
}



/*
 * multiplyRqPolynomials
 * This function returns the coefficients of the R_q multiplied polynomial
 * input of this function is the already calculate sum of logs of the evaluations of the
 * multiplier polynomials. It was chosen this way in order to enable fast computation of
 * the subset-sum, if we use the Gray-code counter mode.
 * Input - array of N uint8_t. the discrete logs of the result polynomial principal
 *				 evaluations on omega**t (when t is odd)
 * Output - the R_q polynomial coefficients
 */
static inline void multiplyRqPolynomials(v16 output[], v8 input[])
{
	// get back from the sum of logs (discrete log) to the multiplication of the principal evaluations
#ifndef v8_signextend_16
	// Use table lookups
	uint8_t *subsetSumOfLogs = (uint8_t*) input;
	for (int i = 0; i < N; ++i) {
		((int16_t*)output)[i] = generatorPowers[subsetSumOfLogs[i]];
	}
#else
	// This is slightly faster because everything is parallelized
	for (int i = 0; i < N/16; ++i) {
		const union cv8 mask15 = {.u8 = {0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
																		 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf}};
		const union cv v2 = {.u16 = {2, 2, 2, 2, 2, 2, 2, 2}};
		v8 x = vector_shuffle(generatorPowersT1, input[i]&mask15.v8);
		// shift right by 4 bits (no instruction with 8-bit words; use 32-bit words)
		v8 t = (v8)v32_shift_r((v32)input[i], 4);
		v8 y = vector_shuffle(generatorPowersT2, t&mask15.v8);
					
		v16 x1 = v8_signextend_16(x);
		v16 y1 = v8_signextend_16(y) + v2.v16;
		output[2*i] = REDUCE(v16_mul(x1, y1));
		x = (v8)v32_shufxor((v32)x, 2); // Swap low half and high half
		y = (v8)v32_shufxor((v32)y, 2); // Swap low half and high half
		x1 = v8_signextend_16(x);
		y1 = v8_signextend_16(y) + v2.v16;
		output[2*i+1] = REDUCE(v16_mul(x1, y1));
	}
#endif

	// change from the point-value evaluation space to the coefficient space
	fft128(output);

	// multiply by powers of OMEGA to get the actual coefficients of the polynomial multiplication
	// apply final modulo. The output is in canonical form, i.e. in the range 0..256
	// the added value CONGRUENCE_ZERO_I was chosen in such a way that the result will always
	// be positive but still it will not cause an overflow of the previous computations
	v16* vecOutputPtr = (v16*)output;

	for (int i = 0; i < N/8; ++i)
	{
		vecOutputPtr[i] = REDUCE(vecOutputPtr[i] * omegaPowers[i]);
	}
}


/*
 * a function to calculate the subset sum of the polynomial exponents, in order to
 * produce the multiplication polynomial results
 * The polynomials in the multiplication are chosen by the boolean array x.
 * One polynomial (index K) is always used in the multiplication.
 * input:		boolean array x of size K
 * output:	subsetSumRq - (size N) the sum of the discrete logs of the of the Rq polynomials
 */
inline void calcSubsetSums (v8 subsetSumRqLogs[], R2poly x)
{
#if SUBSET_TABLES

	union cv8 x8 = {.v8 = (v8)x };

	register v8 Rq_0, Rq_1, Rq_2, Rq_3, Rq_4, Rq_5, Rq_6, Rq_7;

#if PREFETCH
	for (int i=0; i<K/8; i++) {
		// Assuming 64-byte cache line
		__builtin_prefetch(&AA[i][x8.u8[i]][0], 0, 0);
		__builtin_prefetch(&AA[i][x8.u8[i]][4], 0, 0);
	}
#endif

	Rq_0 = AA[0][x8.u8[0]][0];
	Rq_1 = AA[0][x8.u8[0]][1];
	Rq_2 = AA[0][x8.u8[0]][2];
	Rq_3 = AA[0][x8.u8[0]][3];
	Rq_4 = AA[0][x8.u8[0]][4];
	Rq_5 = AA[0][x8.u8[0]][5];
	Rq_6 = AA[0][x8.u8[0]][6];
	Rq_7 = AA[0][x8.u8[0]][7];

	for (int i=1; i<K/8; i++) {
		Rq_0 += AA[i][x8.u8[i]][0];
		Rq_1 += AA[i][x8.u8[i]][1];
		Rq_2 += AA[i][x8.u8[i]][2];
		Rq_3 += AA[i][x8.u8[i]][3];
		Rq_4 += AA[i][x8.u8[i]][4];
		Rq_5 += AA[i][x8.u8[i]][5];
		Rq_6 += AA[i][x8.u8[i]][6];
		Rq_7 += AA[i][x8.u8[i]][7];
	}

	subsetSumRqLogs[0] = Rq_0;
	subsetSumRqLogs[1] = Rq_1;
	subsetSumRqLogs[2] = Rq_2;
	subsetSumRqLogs[3] = Rq_3;
	subsetSumRqLogs[4] = Rq_4;
	subsetSumRqLogs[5] = Rq_5;
	subsetSumRqLogs[6] = Rq_6;
	subsetSumRqLogs[7] = Rq_7;
	
#else
	v8 subsetSumR2Logs [N/32] ALIGN;

	// the polynom whose index is K is always present in the multiplication, regardless of x
	memcpy (subsetSumRqLogs, A[K], N);
	memcpy (subsetSumR2Logs, B[K], N/2);

	v8* vecRqPtr = (v8*)subsetSumRqLogs;

#define likely(x)    __builtin_expect (!!(x), 1)
	
	while (likely(x[0])) {
		int j = __builtin_ffsll(x[0]) - 1;
		x[0] ^= (1ULL<<j);

		// add polynomial j to the subset sum of logs of the R_q polynomial
		v8* vecAptr = (v8*)A[j];
		v8* vecBptr = (v8*)B[j];
		vecRqPtr[0] += vecAptr[0];
		vecRqPtr[1] += vecAptr[1];
		vecRqPtr[2] += vecAptr[2];
		vecRqPtr[3] += vecAptr[3];
		vecRqPtr[4] += vecAptr[4];
		vecRqPtr[5] += vecAptr[5];
		vecRqPtr[6] += vecAptr[6];
		vecRqPtr[7] += vecAptr[7];
	}

#if (K == 128)

	while (likely(x[1])) {
		int j = __builtin_ffsll(x[1]) - 1;
		x[1] ^= (1ULL<<j);

		// add polynomial j to the subset sum of logs of the R_q polynomial
		v8* vecAptr = (v8*)A[j+64];
		v8* vecBptr = (v8*)B[j+64];
		vecRqPtr[0] += vecAptr[0];
		vecRqPtr[1] += vecAptr[1];
		vecRqPtr[2] += vecAptr[2];
		vecRqPtr[3] += vecAptr[3];
		vecRqPtr[4] += vecAptr[4];
		vecRqPtr[5] += vecAptr[5];
		vecRqPtr[6] += vecAptr[6];
		vecRqPtr[7] += vecAptr[7];
	}
#endif

#endif

}

// Extract lsb of representatnt in [-128,128]
// Input is assumed to be in [-127,383]
#ifdef __ARM_NEON__
R2poly getQlsb(const v16 outputQ[])
{
	union u32 polynomQlsb = {.u = {0, 0, 0, 0}};

	for (int i=0; i<16; i+=2) {
		// Put lsb in the sign bit
		// First compare to 128: if lower get lsb, if higher ~lsb
		v16 t1, t2;
		t1 = v16_cmp(outputQ[i  ], V128.v16) ^ v16_shift_l(outputQ[i  ], 15);
		t2 = v16_cmp(outputQ[i+1], V128.v16) ^ v16_shift_l(outputQ[i+1], 15);

		// Then pack 2 vectors of 16-bit integers to 1 vector of 8-bit integers
		// This will overflow, but we only care about the sign bit
		v8 Y = v16_saturate_8(t1,t2);

		// Finally extract sign bit to scalar value (16 bits)
		polynomQlsb.u[i/4] |= v8_extract_msb(Y) << ((i%4)<<3);
	}

	return (v64)polynomQlsb.v;
}
#else
inline R2poly getQlsb(v16 outputQ[])
{
	union u64 polynomQlsb;
	int i;

	// Put lsb in the sign bit
	// First compare to 128: if lower get lsb, if higher ~lsb
	for (i=0; i<N/8; i++)
		outputQ[i] = v16_cmp(outputQ[i], V128.v16) ^ v16_shift_l(outputQ[i], 15);

	// Then pack 2 vectors of 16-bit integers to 1 vector of 8-bit integers
	// This will overflow, but we only care about the sign bit
	v8 Y[N/16];
	for (i=0; i<N/16; i++)
		Y[i] = v16_saturate_8(outputQ[2*i],outputQ[2*i+1]);

	// Finally extract sign bit to scalar value (16 bits)
	for (i=0; i<4; i++)
		polynomQlsb.u[1] = polynomQlsb.u[1]<<16 | v8_extract_msb(Y[7-i]);
	for (i=4; i<8; i++)
		polynomQlsb.u[0] = polynomQlsb.u[0]<<16 | v8_extract_msb(Y[7-i]);

	return polynomQlsb.v;
}
#endif


/*
 * update and modify the running products (resRq, resR2) according to the (one)
 * change in x[] determined by the Gray-code.
 * input: grayCodeCounter - the iteration counter. determines which boolean element in x is changed
 *				x[] - boolean array of size k. determines which polynomials appear in the multiplication
 * modified:	resRq, resR2
 */
void updateXAndProducts (uint64_t* x, v16* resRq, const uint64_t grayCodeCounter)
{
	// determine which element in x is switched
	/* int changeIndex = 0; */
	/* while (!((grayCodeCounter>>changeIndex) & 1))	 { ++changeIndex; } */
	int changeIndex = __builtin_ctz(grayCodeCounter);
	int64_t mask = (1 << changeIndex);

	

	int inv = (*x >> changeIndex) & 1;
	// consider the arrays of values as arrays of SIMD vectors
	v16* vecDRq = (v16*)(D[changeIndex][inv]);

	// update x
	*x ^= mask;

	// update the pointwise product in R_q
	for (int i=0; i<N/8; i++) {
		resRq[i] *= vecDRq[i];
		resRq[i] = REDUCE_FULL_S(resRq[i]);
	}
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
uint64_t BCH128to64 (R2poly packedLSBbits)
{
#if (USE_CLMUL)
	union u64 in = { .v = packedLSBbits };
	// TODO: suboptimal
	union u64 generatorBCHMatrix = {.u = {0xa40137e3da81d585, 0xa40137e3da81d585}};
	R2poly t0 = v64_clmul(packedLSBbits,generatorBCHMatrix.v,0x00);
	R2poly t1 =	v64_clmul(packedLSBbits,generatorBCHMatrix.v,0x11);
	union u64 x = { .v = v64_xor((v64)v32_shufxor((v32)t1,2),t0) };

	return x.u[0] ^ -(in.u[1]&1);
#else
	register uint64_t b1 = packedLSBbits[0];
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

	register uint64_t b2 = packedLSBbits[1];
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
#endif

}

/*
 * apply nIterations of our PRF method in the gray-code counter mode.
 * The results are aggregated somehow into a finalOutput return value,
 * in order to enforce the compiler to do all the instructions and not skip unused ones
 */
uint64_t runGrayCodeMode (uint64_t iv, const int nIterations)
{
	// the subset-sum arrays for both the R_q and the R_2 result polynomials
	v8 subsetSumOfRqLogs [N/16];

	// the calculated coefficients of the R_q result polynomial
	v16 outputQ[N/8];


	// aggregator for the PRF method results. used in order to not enable the compiler
	// to skip some 'unused' code instructions
	uint64_t finalOutput = 0;

	// get the initial subset-sums according to the x boolean array
	union u64 IV = {.v={0, iv}};
	calcSubsetSums (subsetSumOfRqLogs, IV.v);
	// compute the point-wise product
	v16 resRq[N/8];
	int16_t *resRq_scalar = (int16_t*)resRq;
	uint8_t *subsetSumOfRqLogs_scalar = (uint8_t*)subsetSumOfRqLogs;
	for (int i=0; i<N; i++) {
		resRq_scalar[i] = generatorPowers[subsetSumOfRqLogs_scalar[i]];
	}

	// grayCodeCounter is used to determine which boolean element to change in the x array
	uint64_t grayCodeCounter = 0;
	uint64_t x = 0;

	// every output is N/2 bits
	int64_t iter;
	for (iter = 0; iter < nIterations; ++iter)
	{
		// the R_q polynomial multiplication output
		multiplyRqPolynomialsFromPoints (outputQ, resRq);
		// the computation of the round function
		uint64_t res = BCH128to64(getQlsb(outputQ));
		finalOutput ^= res;
		
#if (PRINT_OUTPUT)
		printf ("0x%" PRIx64 " \n", res);
#endif
		// update the gray-code counter and change x and the subset-sums accordingly
		++grayCodeCounter;
		updateXAndProducts (&x, resRq, grayCodeCounter);
	}

	return finalOutput;

}


uint64_t runOutputFeedbackMode (R2poly x, const int nIterations)
{
	if (N != 128 || K != 64)
	{
		fprintf (stderr, "output feedback mode with BCH not implemented for the chosen N,K\n");
		exit(-1);
	}

	// the subset-sum arrays for both the R_q and the R_2 result polynomials
	v8 subsetSumOfRqLogs [N/16];

	// the calculated coefficients of the R_q result polynomial
	v16 outputQ[N/8];

	int64_t iter;
	for (iter = 0; iter < nIterations; ++iter)
	{
		// get the initial subset-sums according to the x boolean array
		calcSubsetSums (subsetSumOfRqLogs, x);
		// the R_q polynomial multiplication output
		multiplyRqPolynomials (outputQ, subsetSumOfRqLogs);
		
		x = getQlsb(outputQ);
		x[0] = BCH128to64(x);

#if (PRINT_OUTPUT)
		printf ("0x%" PRIx64 " \n", (uint64_t)(x[0]));
#endif
	}

	return x[0];
}

int main ()
{

	// seed the key - initialize the polynomial tables
	initializeA();
	initializeB_C();

	// compute the pointwise values of the polynomial
	// this require generatorPowers to be initialized
	initializeD();

	// intialize precomputed AA
	initializeXX();

	printf ("SPRING-BCH: N = %d	 K = %d\n\n", N, K);


	uint64_t iv;
	// iv = 0xaaaaaaaabbbbbbbb;
	iv = 0;


#ifdef rdtsc
	uint64_t tsc = rdtsc();
#else
	clock_t begin, end;
	begin = clock();
#endif

	uint64_t finalOutput;
	finalOutput = runGrayCodeMode (iv, N_ITERATIONS);

#ifdef rdtsc
	tsc = rdtsc() - tsc;
#else
	end = clock();
#endif

	printf ("COUNTER MODE:\n");
	printf ("output = %" PRIx64 " \n", finalOutput);
#ifdef rdtsc
	printf ("%f c/B\n", 8.*tsc/(1.*N_ITERATIONS*(N/2)));
#else
	double dt = (double) (end - begin) / CLOCKS_PER_SEC;
	printf ("%f MB/s (time = %f)\n", ((float)N_ITERATIONS*(N/2)/8000000)/dt, dt);
#endif
	printf ("\n");

#if K == 64

	union u64 x;
	//x.u[0] = x.u[1] = 0xaaaaaaaabbbbbbbb;
	x.u[0] = x.u[1] = 0;

#ifdef rdtsc
	tsc = rdtsc();
#else
	begin = clock();
#endif

	finalOutput = runOutputFeedbackMode (x.v, N_ITERATIONS);

#ifdef rdtsc
	tsc = rdtsc() - tsc;
#else
	end = clock();
	dt = (double) (end - begin) / CLOCKS_PER_SEC;
#endif

	printf ("OFB MODE: \n");
	printf ("output = %" PRIx64 " \n", finalOutput);
#ifdef rdtsc
	printf ("%f c/B\n", 8.*tsc/(1.*N_ITERATIONS*(N/2)));
#else
	printf ("%f MB/s (time = %f)\n", ((float)N_ITERATIONS*(N/2)/8000000)/dt, dt);
#endif
	printf ("\n");

#endif

	return 0;
}

// Local Variables: 
// c-file-style: "linux"
// c-basic-offset: 2
// tab-width: 2
// indent-tabs-mode: t
// End: 
