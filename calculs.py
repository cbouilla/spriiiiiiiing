from fractions import Fraction
from math import log

def fact(n):
	f = 1
	for i in range(1, n+1):
	     f *= i
	return f

def binomial(n, k):
	return fact(n) // (fact(k) * fact(n-k))

def term_k(n, i, p=Fraction(1, 257)):
	return binomial(n, i) * pow(p, i) * pow(1-p, n-i)

def f(n, k, p=Fraction(1, 257)):
    """proba que >= k coeffs soient nuls sur n. Un coeff est nul avec proba p"""
    return 1 - sum([term_k(n, i, p) for i in range(k+1)])


#for i in range(32):
#	x = f(128, i)
#	print("proba d'avoir > {0} bit foireux parmi 128 : {1}".format(i, log(x.numerator, 2) - log(x.denominator, 2)))


w = 0
for i in range(1, 9):
	w += (8-i) * term_k(8, i)

print("SSE : {0}%".format( float(w)/8*100)) # average number of wasted bits


w = 0
for i in range(1, 17):
	w += (16-i) * term_k(16, i)

print("AVX2 : {0}%".format( float(w)/16*100)) # average number of wasted bits

w = 0
for i in range(1, 33):
	w += (32-i) * term_k(32, i)

print("AVX512 : {0}%".format( float(w)/32*100)) # average number of wasted bits

