
p = 2
Z_p = IntegerModRing(p)
Z_pX, lamb = Z_p["λ"].objgen()

L = 2
currMatrixElems = [[0 for x in range(0, L)] for y in range(0, L)]
rolledOver = False

numFactorable = 0

while not rolledOver:
    currMatrix = matrix(Z_p, currMatrixElems)
    factors = list(currMatrix.charpoly().factor())
    
    if (len(factors) > 1):
        numFactorable += 1
    else: #Check for repeated factors
        for f in factors:
            if (f[1] > 1):
                numFactorable += 1
                break
    
    #Increment currMatrix
    rolledOver = True
    for x in range(0, L):
        for y in range(0, L):
            if (currMatrixElems[x][y] == p-1):
                currMatrixElems[x][y] = 0
            else:
                currMatrixElems[x][y] += 1
                rolledOver = False
                break
        if not rolledOver:
            break
        
print(f"Number of matrices mod {p} with factorable characteristic polynomials: {numFactorable}")
print(f"Predicted bounds: [{3*p^3 - 4*p^2 + 3*p}, {p^4 - p^3 + p^2 + p}]")
print(f"Predicted count using average: {(p^4 + 2*p^3 - p^2)/2}")
print(f"Difference between real count and prediction: {numFactorable - (p^4 + 2*p^3 - 3*p^2 + 4*p)/2}")

#My original formula was off by p(p-2) ... (3,15,35,99,143,255,323...)

'''
Collected data:

modulus used             : 2  3  5   7    11   13    17    19    23     29
# factorable chara polys : 14 63 425 1519 8591 16393 46529 71839 151823 377609
formula estimation       : 14 60 410 1484 8492 16250 46274 71516 151340 376826
'''