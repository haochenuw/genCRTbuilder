load('utils.py')

class CRTbuilder:
    
    def __init__(self, m,t):
        """
        cosetrepMap -- maps integer to its canonical coset in (Zm)*/<t>.
        """
        self.m = m
        self.t = t
        self.numslots = get_numslots(self.m, self.t)
        self.rootTable = compute_rootTable(self.m, self.numslots, self.t)
        self.cosetreps = computecosetRep_new(self.m,self.t)
        self.n = euler_phi(self.m)
        
        self.baseField = self.rootTable[0].parent()
        self.invrootTable = [a**(-1) for a in self.rootTable]
        self.logn = round(log(self.n,2))
        
        self.primeField = self.baseField.prime_subfield()
        
        self.frob = self.baseField.frobenius_endomorphism()

    def decompose(self, data):
        """
        return a dictionary
        """
        # this part data should have length n
        result = []
        if len(data) != self.n:
            raise ValueError
        try:
            for i in range(self.n):
                data[i] = self.baseField(data[i])
        except:
            raise ValueError("plaintext do not lie in appropriate field.")
            
            
        output = negacyclic_ntt(data, self.invrootTable)
        for j in range(self.n):
            power = 2*j+1
            if self.cosetreps[power][0] == power:
                result.append(output[j])
        return result
        
    def compose(self, data):
        # ensure that data has correct length
        d = len(data)
        if d != self.numslots:
            raise ValueError
        
        classes = self.cosetreps
        # expand back
        expanded = []
        for j in range(self.n):
            power =  2*j + 1
            reducedpower = classes[power][0]
            j1 = (reducedpower -1) // 2
            try:
                expanded.append((self.frob**classes[power][1])(data[j1]))
            except:
                raise ValueError(power,classes[power])# todo: add frobenius,
        result =  inv_negacyclic_ntt(expanded, self.rootTable)
        assert len(result) == self.n
        for i in range(self.n):
            if result[i] not in self.primeField:
                raise ValueError("sth wrong")         
        return result
            
            
        
    
        
        