from random import randrange, seed

seed(int(12550821))

N = 49
Z_N = IntegerModRing(N)

def vec_span(modulus, listOfVects):
    '''Because Sage is grumpy and won't compute spans for moduli
       which aren't prime...
       
       The function assumes all vectors in listOfVects are the same size.'''
    integerRing = IntegerModRing(modulus)
    coeffs = [0 for x in range(0, len(listOfVects))]
    linearCombos = []
    
    #Creating every possible linear combination of vectors in listOfVects
    while (coeffs[-1] < modulus):
        tempVect = vector(integerRing, [0 for x in range(0, len(listOfVects[0]))])
        
        #Create linear combination
        for c, v in enumerate(listOfVects):
            tempVect += coeffs[c]*v
            
        if (tempVect not in linearCombos):
            linearCombos.append(tempVect)
            
        #Increment coeffs
        for i in range(0, len(coeffs)):
            coeffs[i] += 1
            if ((coeffs[i] >= modulus) and (i != len(coeffs)-1)):
                coeffs[i] = 0
            else:
                break
                
    return linearCombos


#Generate random matrix
L = 5
tempMatrix = []
for row in range(0, L):
    tempMatrix.append([randrange(0, N) for x in range(0, L)])
A = matrix(Z_N, tempMatrix)
print(A)

#Create our vector space
'''
basisVectors = []
for i in range(0, L):
    tempVect = [0 for x in range(0, L)]
    tempVect[i] = 1
    basisVectors.append(vector(Z_N, tempVect))
basisSpace = vec_span(N, basisVectors)
print(basisSpace)
'''