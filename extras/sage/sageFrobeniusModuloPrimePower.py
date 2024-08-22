
from random import randrange, seed
seed(int(12550821))

p = 7
k = 2
Z_pk = IntegerModRing(p^k)

L = 5 #Size of matrices

#The "Frobenius Normal Form" we're using
F = matrix(Z_pk, [[3, 1, 0, 0, 0],
                  [2, 0, 1, 0, 0],
                  [1, 0, 0, 0, 0],
                  [0, 0, 0, 6, 1],
                  [0, 0, 0, 33, 0]])

#Create a random change of basis matrix
B = None
while True:
    tempMatrix = [[randrange(0, p^k) for x in range(0, L)] for y in range(0, L)]
    B = matrix(Z_pk, tempMatrix)
    if (B.charpoly()[0] % p != 0):
        break
        
A = B*F*B^(-1)
print("A's characteristic polynomial:", A.charpoly())

print("F = ", F, sep="\n", end="\n\n")
print("BFB^(-1) = ", A, sep="\n", end="\n\n")

#Creating vectors to see if A has an "easy" cyclic decomposition
# using the vectors that work for F
v_1 = B*vector(Z_pk, [0, 0, 1, 0, 0])*B^(-1)
v_2 = B*vector(Z_pk, [0, 0, 0, 0, 1])*B^(-1)

print("v_1's orbit:")
for i in range(0, L):
    print((A^i)*v_1)
    
print("v_2's orbit:")
x = []
for i in range(0, L):
    print((A^i)*v_2)
    x.append((A^i)*v_2)
    
#Verifying that vectors are, indeed, linearly dependent
print("")
coeffs = [0 for i in range(0, L)]
specialCoeff = 0
rolledOver = False
while not rolledOver:
    coeffs[specialCoeff] = 1
    tempVector = vector(Z_pk, [0 for i in range(0, 5)])
    for i in range(0, L):
        tempVector += x[i]*coeffs[i]
        
    if (tempVector == 0):
        print("With specialCoeff =", specialCoeff, ":", coeffs)
        specialCoeff += 1
        
    #Increment coeffs
    rolledOver = True
    if (specialCoeff < L):
        for i in range(0, L):
            if (i == specialCoeff):
                continue
            coeffs[i] += 1
            if (coeffs[i] >= p^k):
                coeffs[i] = 0
            else:
                rolledOver = False
                break
                
    if ((rolledOver) and (specialCoeff < L)):
        specialCoeff += 1
        roleldOver = False
