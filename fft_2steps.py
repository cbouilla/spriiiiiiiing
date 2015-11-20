import random

# la taille du corps
q = 257
log_n = 4
n = 1 << log_n

class ShortRange:
    a = -128 # min
    b = 383  # max

    def __add__(self, other):
        result = ShortRange()
        result.a = self.a + other.a
        result.b = self.b + other.b
        return result

    def __sub__(self, other):
        result = ShortRange()
        result.a = self.a - other.b
        result.b = self.b - other.a
        return result

    def __mul__(self, k):
        result = ShortRange()
        result.a = k * self.a
        result.b = k * self.b
        return result

    def overflow(self):
        return (self.a < -32768) or (self.b > 32767)


root = {256: 3, 
        128: 42,
        64: -35,
        32: -60,
        16: 2, 
        8:  4, 
        4: 16, 
        2: 256}
cst_i = 16 # cst_i**2 == -1

def slow_dft(a, _n, omega):
    b = [0 for i in range(_n)]
    for i in range(_n):
        x = pow(omega, i, q)
        for j in range(_n):
            b[i] = (b[i] + a[j] * pow(x, j, q)) % q
    return b

# renvoie le nombre obtenu en écrivant les bits de x en sens inverse
def bit_reverse(x, n):
    if n == 2:
        return x
    if n == 4:
        return ((x & 0x01) << 1) | ((x & 0x02) >> 1)
    if n == 8:
        return ((x & 0x01) << 2) | (x & 0x02) | ((x & 0x04) >> 2)
    if n == 16:
        x = ((x & 0x05) << 1) | ((x & 0x0A) >> 1)
        x = ((x & 0x03) << 2) | ((x & 0x0C) >> 2)
        return x
    assert False



# decimation-in-frequency
# input in normal order
# output in bit-reversed-order
def dif_butterfly(a, i, j, omega):
    global adds, mults
    #print("  DIF-butterfly {0} <---> {1}  [omega={2}]".format(i, j, omega))
    u = (a[i] + a[j]) % q
    v = ((a[i] - a[j]) *  omega) % q
    a[i] = u
    a[j] = v
    adds += 2
    if omega != 1 and omega != 256:
        mults += 1


def dif(a, k, start, omega, stride=1):
    """
    omega est une racine k-ème de l'unité
    """
    if k == 1:
        return

    #print("size-{0} DIF-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
    for i in range(k//2):
        dif_butterfly(a, start + stride*i, start + stride*(i + k//2), pow(omega, i))
    omega_prime = pow(omega, 2, q)
    dif(a, k//2, start, omega_prime, stride)
    dif(a, k//2, start + stride*(k//2), omega_prime, stride)


# decimation-in-time
# input in bit-reversed order
# output in normal order
def dit_butterfly(a, i, j, omega):
    global adds, mults

    #print("  DIT-butterfly {0} <---> {1}  [omega={2}]".format(i, j, omega))
    u = a[i]
    v = (a[j] *  omega) % q
    a[i] = (u + v) % q
    a[j] = (u - v) % q
    adds += 2
    if omega != 1 and omega != 256:
        mults += 1


def dit(a, k, start, omega, stride=1):
    """
    omega est une racine k-ème de l'unité
    """
    if k == 1:
        return

    #print("size-{0} DIT-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
    omega_prime = pow(omega, 2, q)
    dit(a, k//2, start, omega_prime, stride)
    dit(a, k//2, start + stride*(k//2), omega_prime, stride)
    for i in range(k//2):
        dit_butterfly(a, start + stride*i, start + stride*(i + k//2), pow(omega, i))


def center(x):
    x = x % q
    if x > (q // 2):
        x -= q
    if x < -(q // 2):
        x += q
    return x


adds = 0
mults = 0

a = 8
b = 16
N = a*b
omega = root[N]
A = [ center(5*i ^ i*17 ^ 42) for i in range(N)]


print("input data : ")
for i in range(a):
    for j in range(b):
        print( '{0:4d}, '.format(A[b*i + j]), end="" )
    print()


B = slow_dft(A, N, omega)

# matrix is : a (rows) * b (columns)

# b * size-a FFT (all columns in parallel)
for i in range(b):
    dif(A, a, i, pow(omega, b, q), stride=b)

print("after step 1 : ")
for i in range(a):
    for j in range(b):
        print( '{0:4d}, '.format(center(A[b*i + j])), end="" )
    print()


# twiddle
T = []
for i in range(a):
    for j in range(b):
        k = bit_reverse(i, a) * j
        A[b*i + j] *= pow(omega, k, q)
        T.append(pow(omega, k, q));

print("after twiddle : ")
for i in range(a):
    for j in range(b):
        print( '{0:4d}, '.format(center(A[b*i + j])), end="" )
    print()


# transpose
C = []
for i in range(a):
    for j in range(b):
        C.append( A[ bit_reverse(i, a) + a*bit_reverse(j, b) ] )
A = C

print("after transpose : ")
for i in range(a):
    for j in range(b):
        print( '{0:4d}, '.format(center(A[b*i + j])), end="" )
    print()

# matrix is now : b (rows) * a (columns)
# a * size-b FFT  (all columns in parallel)
for i in range(a):
    dit(A, b, i, pow(omega, a, q), stride=a)

print("after final step : ")
for i in range(a):
    for j in range(b):
        print( '{0:4d}, '.format(center(A[b*i + j])), end="" )
    print()


for i in range(N):
    assert A[i] == B[i]