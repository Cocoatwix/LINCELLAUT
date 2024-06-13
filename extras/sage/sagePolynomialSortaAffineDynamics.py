MODULUS = 41
Z_N = IntegerModRing(MODULUS)
polysOverZ_N, lamb = Z_N["λ"].objgen()

#The minimal polynomial defining our quotient ring
# polysOverZ_N / MINPOLY
MINPOLY = 8*lamb - 8

#Formula for the order of Σ_i x^i over polysOverZ_N / aX - c
#print( ( Z_N(c) * (Z_N(a)^(-1)) ).multiplicative_order())

POWERSPACING = 1 #How much the exponents of the different powers of λ we're adding differ by
runningSum   = 1 #Holds the polynomial sum we're interested in
currentPower = 1 #Holds the current power we're adding to runningSum

order = 0 #Keeps track of how many iterations it takes for Σ x^{POWERSPACING*i} to iterate back to 1

started = False #For simulating a do-while loop

#Loop until our runningSum loops back to 1
print("+ | Σ")
while ((runningSum != 1) or (currentPower != 1) or (not started)):
    started = True
    
    order += 1
    
    currentPower *= lamb^POWERSPACING
    currentPower %= MINPOLY
    runningSum += currentPower
    runningSum %= MINPOLY
    print(currentPower, runningSum, sep=" ; ")
    
print("Order:", order)
print("---------")
