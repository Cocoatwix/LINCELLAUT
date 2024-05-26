Z_7 = IntegerModRing(7)
ZP_7, lamb = Z_7["λ"].objgen() #Polynomial ring over Z_7
# f = 2*lamb^2 + 1

R = matrix(Z_7, [[1, 2, 0, 6], 
                 [0, 1, 2, 3], 
                 [0, 0, 1, 1], 
                 [0, 0, 0, 1]])
A = matrix(Z_7, [[3, 1, 0, 0], 
                 [2, 0, 1, 0],
                 [5, 0, 0, 1],
                 [5, 0, 0, 0]])
B = (R^(-1))*A*R
I = matrix(Z_7, [[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0],
                 [0, 0, 0, 1]])

#Finding the multiplicative order of B
i = 1
while (B^i != I):
    i += 1
print(i)

#factor(B.minpoly("λ"))