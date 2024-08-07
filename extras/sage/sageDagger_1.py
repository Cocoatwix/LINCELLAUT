
#https://doc.sagemath.org/html/en/reference

#Setting up the modulus
p = 7
k = 2
N = p^k
Z_p = IntegerModRing(p)
Z_N = IntegerModRing(N)

def dagger_1(A, C, n):
    '''Computes dagger_1^n(A, C).'''
    tempSum = matrix(Z_N, len(A.rows()), len(A.rows()[0]))
    
    for i in range(0, n):
        tempSum += (A^i)*(C)*(A^(n-i-1))
        
    return tempSum
	
	
'''
To get the dagger matrix, I see where matrices with a 1 in one component and 0s everywhere else
get sent (the "basis vectors" for the matrices), then I convert those resulting matrices into column 
vectors by stacking the column vectors of the matrices together, leftmost vector at the top, 
rightmost at the bottom.

For example, to get the first column vector for the dagger matrix where
A = Fibonacci matrix, C = 1 in top-left corner and zero everywhere else,
and n = 2, I compute AC + CA (mod 5), which gives me
[2 1]
[1 0].
The column vector for the dagger matrix is then
[2]
[1]
[1]
[0].
'''
def dagger_1_matrix(A, n):
	'''Create the dagger matrix representation for dagger_1^n(A).'''
	rowNumA = len(A.rows())
	colNumA = len(A.rows()[0])

	daggerMatrixElements = [[0 for col in range(0, colNumA*rowNumA)] for row in range(0, rowNumA*colNumA)]
	tempComponents = [[0 for col in range(0, colNumA)] for row in range(0, rowNumA)]

	#Compute dagger matrices for each "basis matrix"
	for i in range(0, rowNumA*colNumA):
		if (i != 0):
			tempComponents[(i-1) % rowNumA][(i-1) // colNumA] = 0
		
		tempComponents[i % rowNumA][i // colNumA] = 1
		tempComponent = matrix(Z_N, tempComponents)
		tempDagger = dagger_1(A, tempComponent, n)
		
		#Now, convert the computed dagger matrix into a column vector
		for col in range(0, colNumA):
			for row in range(0, rowNumA):
				daggerMatrixElements[row + col*rowNumA][i] = tempDagger[row][col]
				
	return matrix(Z_N, daggerMatrixElements)

matrixOfInterest = [[1, 1],
                    [0, 1]]
tau = 0 #Transient length of A
A = matrix(Z_N, matrixOfInterest)
A_p = matrix(Z_p, matrixOfInterest)
C = matrix(Z_N, [[3, 3], [1, 4]])
n = 7

print("Minimal polynomial of A:", A_p.minpoly(), end="\n\n")
print("A^{" + str(n) + "} - A^{tau} mod " + str(N) + ":", A^n - A^tau, sep="\n")
print("\n")
print("(A + " + str(p) + "^{" + str(k-1) + "}C)^" + str(n) + " - A^{tau} mod " + str(N) + ":", (A + (p^(k-1))*C)^n - A^tau, sep="\n")
print("\n")
print("Dagger \"basis\" matrix:", (p^(k-1))*dagger_1_matrix(A, n), sep="\n")
