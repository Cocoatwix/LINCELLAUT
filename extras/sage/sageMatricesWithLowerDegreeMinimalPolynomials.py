from random import randrange, seed

seed(int(12550821))

p = 7
k = 2
N = p^k
Z_p = IntegerModRing(p)
Z_N = IntegerModRing(N)

#Generate a random matrix where the minimal polynomial is
# of lesser degree than the characteristic.
L = 5
while True:
    tempMatrix = []
    for row in range(0, L):
        tempMatrix.append([randrange(0, N) for x in range(0, L)])
    A_p = matrix(Z_p, tempMatrix)
    A   = matrix(Z_N, tempMatrix)
    if (A_p.minpoly().degree() < L):
        break
 
#Print the matrix without surrounding brackets
for row in range(0, L):
    for col in range(0, L):
        print(A[row][col], end=" ")
    print("")
