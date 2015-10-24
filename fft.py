import random

# la taille du corps
q = 257
log_n = 4
n = 1 << log_n

root = {16: 2, 
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
def bit_reverse(x):
    if n == 4:
        return ((x & 0x01) << 1) | ((x & 0x02) >> 1)
    if n == 8:
        return ((x & 0x01) << 2) | (x & 0x02) | ((x & 0x04) >> 2)
    if n == 16:
        x = ((x & 0x05) << 1) | ((x & 0x0A) >> 1)
        x = ((x & 0x03) << 2) | ((x & 0x0C) >> 2)
        return x
    assert False

def test_bit_reverse():
    assert bit_reverse(0) == 0
    for i in range(log_n):
        assert bit_reverse(1 << i) == 1 << (log_n-1-i)

    for i in range(n):
        assert bit_reverse(bit_reverse(i)) == i

    for i in range(n):
        for j in range(n):
            assert bit_reverse(i ^ j) == bit_reverse(i) ^ bit_reverse(j)

# renvoie le nombre obtenu en écrivant les bits de x en sens inverse
def bit_reverse_radix4(x):
    if n == 16:
        return ((x & 0x03) << 2) | ((x & 0xc) >> 2)
    assert False

def test_bit_reverse_radix4():
    assert bit_reverse_radix4(0) == 0
    assert bit_reverse_radix4(1) == 4
    assert bit_reverse_radix4(2) == 8

    for i in range(n):
        assert bit_reverse_radix4(bit_reverse_radix4(i)) == i

    for i in range(n):
        for j in range(n):
            assert bit_reverse_radix4(i ^ j) == bit_reverse_radix4(i) ^ bit_reverse_radix4(j)


# decimation-in-frequency
# input in normal order
# output in bit-reversed-order
def dif_butterfly(a, i, j, omega):
    global adds, mults
    print("  DIF-butterfly {0} <---> {1}  [omega={2}]".format(i, j, omega))
    u = (a[i] + a[j]) % q
    v = ((a[i] - a[j]) *  omega) % q
    a[i] = u
    a[j] = v
    adds += 2
    if omega != 1 and omega != 256:
        mults += 1

# computes a DFT-4
def dif_butterfly_radix4(a, i, j, k, l, omega):
    global adds, mults
    print("  DIF-butterfly-4 {0} [omega={1}]".format([i, j, k, l], omega))
    
    
    #DIF-butterfly 0 <---> 2  [omega=1]
    u = a[i] + a[k]
    w = a[i] - a[k]

    #DIF-butterfly 1 <---> 3  [omega=16]
    v =  a[j] + a[l]
    x = (a[j] - a[l]) * cst_i

    # j <---> k à cause du bit reversal pour avoir la sortie dans le bon ordre
    #DIF-butterfly 0 <---> 1  [omega=1]
    a[i] =  u + v
    a[k] = (u - v) * omega**2 
    
    #DIF-butterfly 2 <---> 3  [omega=1]
    a[j] = (w + x) * omega
    a[l] = (w - x) * omega**3

    a[i] %= 257
    a[j] %= 257
    a[k] %= 257
    a[l] %= 257

    adds += 8
    mults += 1
    if omega != 1 and omega != 256:
        mults += 1
    if pow(omega, 2, q) != 1 and pow(omega, 2, q) != 256:
        mults += 1
    if pow(omega, 3, q) != 1 and pow(omega, 3, q) != 256:
        mults += 1

def test_radix4_butterfly():
    a = [random.randrange(q) for i in range(4)]
    b = slow_dft(a, 4, cst_i)
    dif_butterfly_radix4(a, 0, 1, 2, 3, 1)
    for i in range(4):
        assert a[i] == b[i]


def dif(a, k, start, omega):
    """
    omega est une racine k-ème de l'unité
    """
    if k == 1:
        return

    #print("size-{0} DIF-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
    for i in range(k//2):
        dif_butterfly(a, start + i, start + i + k//2, pow(omega, i))
    omega_prime = pow(omega, 2, q)
    dif(a, k//2, start, omega_prime)
    dif(a, k//2, start + k//2, omega_prime)

def test_dif():
    omega = root[n]
    a = [random.randrange(q) for i in range(n)]
    b = slow_dft(a, n, omega)
    dif(a, n, 0, omega)
    for i in range(n):
        assert a[bit_reverse(i)] == b[i]


def dif_radix4(a, k, start, omega):
    """
    omega est une racine k-ème de l'unité
    """
    if k == 1:
        return
    print("size-{0} DIF-radix4-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
    for i in range(k//4):
        dif_butterfly_radix4(a, start + i, start + i + k//4, start + i + k//2, start + i + 3*k//4, pow(omega, i))
    omega_prime = pow(omega, 4, q)
    for i in range(k//4):
        dif_radix4(a, k//4, start + i*k//4, omega_prime)

def test_dif_radix4():
    omega = root[n]
    a = [random.randrange(q) for i in range(n)]
    b = slow_dft(a, n, omega)
    dif_radix4(a, n, 0, omega)
    print(a)
    print(b)
    for i in range(n):
        assert a[bit_reverse_radix4(i)] == b[i]


# decimation-in-time
# input in bit-reversed order
# output in normal order
def dit_butterfly(a, i, j, omega):
    global adds, mults

    print("  DIT-butterfly {0} <---> {1}  [omega={2}]".format(i, j, omega))
    u = a[i]
    v = (a[j] *  omega) % q
    a[i] = (u + v) % q
    a[j] = (u - v) % q
    adds += 2
    if omega != 1 and omega != 256:
        mults += 1


def dit(a, k, start, omega):
    """
    omega est une racine k-ème de l'unité
    """
    if k == 1:
        return

    #print("size-{0} DIT-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
    omega_prime = pow(omega, 2, q)
    dit(a, k//2, start, omega_prime)
    dit(a, k//2, start + k//2, omega_prime)
    for i in range(k//2):
        dit_butterfly(a, start + i, start + i + k//2, pow(omega, i))

def test_dit():
    omega = root[n] # très pratique
    a = [random.randrange(q) for i in range(n)]
    b = slow_dft([a[bit_reverse(i)] for i in range(n)], n, omega)
    dit(a, n, 0, omega)
#    print(a)
#    print(b)
    for i in range(n):
        assert a[i] == b[i]
    

def split_dit_butterfly(a, i, j, k, l, omega):
    global adds, mults
    print("  split-DIT-butterfly E={0}, {1} / O = {2}, {3} [omega={2}]".format(i, j, k, l, omega))
    e0 = a[i]
    e1 = a[j]
    o0 = a[k] * omega
    o1 = a[l] * pow(omega, 3)
    A = o0 + o1
    B = (o0 - o1) * cst_i
    a[i] = (e0 + A) % q
    a[k] = (e0 - A) % q
    a[j] = (e1 + B) % q
    a[l] = (e1 - B) % q
    adds += 6
    mults += 1
    if omega != 1 and omega != 256:
        mults += 1
    if pow(omega,3, q) != 1 and pow(omega, 3, q) != 256:
        mults += 1

def test_splitradix_butterfly():
    a = [random.randrange(q) for i in range(4)]
    b = slow_dft(a, 4, cst_i)
    x = a[0]
    u = a[1]
    y = a[2]
    v = a[3]
    a = [x + y, x-y, u, v]
    split_dit_butterfly(a, 0, 1, 2, 3, 1)
    print(a)
    print(b)
    for i in range(4):
        assert a[i] == b[i]

# output in the right order
# input in the "split-radix" revbin order
def splitradix_dit(a, k, start, omega):
    """
    omega est une racine k-ème de l'unité
    """
    if k == 1:
        return
    if k == 2:
        print("size-{0} split-radix DIT-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
        u = a[start]
        v = a[start+1]
        a[start  ] = (u + v) % q
        a[start+1] = (u - v) % q
        return

    print("size-{0} split-radix DIT-FFT on [{1}:{2}] with omega={3}".format(k, start, start+k, omega))
    
    splitradix_dit(a, k//2, start,          pow(omega, 2, q))
    splitradix_dit(a, k//4, start +   k//2, pow(omega, 4, q))
    splitradix_dit(a, k//4, start + 3*k//4, pow(omega, 4, q))

    for i in range(k//4):
        split_dit_butterfly(a, start+i, start+i + k//4, start+i+k//2, start+i+3*k//4, pow(omega, i))

def test_splitradix_dit():
    omega = root[n] # très pratique
    a = [random.randrange(q) for i in range(n)]
    b = slow_dft(a, n, omega)
    # put this mess in the right order
    if n == 4:
        a = [a[0], a[2], a[1], a[3]]
    elif n == 8:
        a = [a[0], a[4], a[2], a[6], a[1], a[5], a[3], a[7]]
    elif n == 16:
        a = [a[0], a[8], a[4], a[12], a[2], a[10], a[6], a[14], a[1], a[9], a[5], a[13], a[3], a[11], a[7], a[15]]
    else:
        raise ValueError("unsupported size")
    splitradix_dit(a, n, 0, omega)
    #print(a)
    #print(b)
    for i in range(n):
        assert a[i] == b[i]
    


#test_bit_reverse()
#test_bit_reverse_radix4()

adds = 0
mults = 0
#test_radix4_butterfly()
#test_splitradix_butterfly()
test_splitradix_dit()
print("split-radix = {0} additions, {1} multiplications".format(adds, mults))   
print("------------------------------") 

adds = 0
mults = 0
test_dif()
print("radix-2 = {0} additions, {1} multiplications".format(adds, mults))   
print("------------------------------") 

if n == 16:
    adds = 0
    mults = 0
    test_dif_radix4()
    print("radix-4 = {0} additions, {1} multiplications".format(adds, mults))   

