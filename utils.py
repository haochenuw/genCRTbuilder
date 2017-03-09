

"""
reversing the bits and reconstruct the integer
"""
def reverse_bits(x, Maxbits):
    res = x.bits() 
    res += [0 for _ in range(Maxbits - len(res))]
    res.reverse()
    return sum([res[i]*2**i for i in range(Maxbits)])



def get_numslots(m,t):
    return euler_phi(m) // Zmod(m)(t).multiplicative_order()


def compute_rootTable(m,numslots,t):
    if not is_prime(t):
        raise ValueError
    n = euler_phi(m)
    # check that n is a power of 2
    ext_deg = n // numslots
    FF = FiniteField(t**ext_deg)
    psi = FF.multiplicative_generator()**((FF.cardinality() - 1) // m)
    assert psi.multiplicative_order() == m
    return [psi**reverse_bits(Integer(j),log(n,2)) for j in range(n)]


"""
compute coset representatives for Zm*/(t), plus the number of hops
"""
def computecosetRep(m,t):
    if gcd(m,t) != 1:
        raise ValueError
    classes = [0 for _ in range(m)]
    for i in range(m):
        if gcd(i,m) == 1:
            classes[i] = i
    for i in range(m):
        if classes[i] ==0: continue

        if classes[i] < i:
            classes[i] = classes[classes[i]]
            continue

        j = (i*t) % m
        while classes[j] != i:
            classes[classes[j]] = i
            j = (j*t) % m
    return classes


def computecosetRep_new(m,t):
    if gcd(m,t) != 1:
        raise ValueError
    classes = []
    for i in range(m):
        if gcd(i,m) == 1:
            classes.append([i,0])
        else:
            classes.append([0,0])
    for i in range(m):
        if classes[i][0] ==0: 
            continue
        if classes[i][0] < i:
            continue
        j = (i*t) % m
        hop = 1
        while classes[j][0] != i:
            classes[j][0] = i
            classes[j][1] = hop
            j = (j*t) % m
            hop +=1
    return classes



def negacyclic_ntt(data, rootTable): 
    n = len(data)
    tt = n
    FF = rootTable[0].parent()
    for k in range(log(n,2)):
        mm = 2**k
        tt /= 2
        for i in range(mm):
            j1 = 2*i*tt
            j2 = j1+tt-1
            S = rootTable[mm+i]
            for j in range(j1,j2+1):
                U = FF(data[j])
                V = FF(data[j+tt]*S)
                data[j] = U + V
                data[j+tt] = U - V
    return data
                


def inv_negacyclic_ntt(data, invrootTable):
    tt = 1
    n = len(data)
    logn = round(log(n,2))
    for k in range(logn):
        mm = 2**(logn-k)
        j1=  0
        h = mm /2
        for i in range(h):
            j2 = j1+tt-1
            S = invrootTable[h+i]
            for j in range(j1, j2+1):
                U,V = data[j], data[j+tt]
                data[j] = U+V
                data[j+tt] = (U-V)*S
            j1 += 2*tt
        tt *= 2
    baseField = invrootTable[0].parent()
    ninv = baseField(n)**(-1)
    for j in range(n):
        data[j] *= ninv
    return data

