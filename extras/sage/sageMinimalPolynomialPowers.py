from random import randrange, seed
#seed(int(12550821)) #Integers by default are not regular Python integers in Sage

MODULUS = 23
Z_N = IntegerModRing(MODULUS)
Z_NX, lamb = Z_N["λ"].objgen()

#print((lamb^(MODULUS^3) - lamb).factor())

definingMinPolyFactor = lamb^3 + 7*lamb^2 + 15*lamb + 7
definingMinPolyMult   = 3
definingMinPoly       = definingMinPolyFactor^definingMinPolyMult
print(f"Modulus: {MODULUS}")
print(f"Defining minimal polynomial: ({definingMinPolyFactor})^{definingMinPolyMult} = {definingMinPoly}")

definingCoeffs  = [] #Holds proper min poly coeffs for creating companion matrix
for i in range(0, definingMinPoly.degree()):
    definingCoeffs.append(MODULUS - definingMinPoly[i])

#Creating our companion matrix
#As well, creating our random det 1 matrix to transform the
# companion matrix into something more arbitrary.
companElems = []
detElems    = []

for row in range(0, definingMinPoly.degree()):
    companElems.append([definingCoeffs[-row-1]])
    companElems[-1] += [0]*row
    
    detElems.append([0]*row)
    detElems[-1] += [1]
    detElems[-1] += [randrange(0, MODULUS) for e in range(row+1, definingMinPoly.degree())]
    
    if (definingMinPoly.degree() - row - 2 >= 0):
        companElems[-1] += [1]
        companElems[-1] += [0]*(definingMinPoly.degree()-row-2)
        
companionMatrix = matrix(Z_N, companElems)
detMatrix       = matrix(Z_N, detElems)

transformedMatrix = (detMatrix^(-1))*companionMatrix*(detMatrix)
print(f"Transformed matrix:\n{transformedMatrix}")

#Now, look at the sequence of kernels created by
# increasing the multiplicity of the irreducible factor...
for i in range(1, definingMinPolyMult+1):
    print(f"{i}:\n", kernel((definingMinPolyFactor^i)(transformedMatrix)), sep="")
