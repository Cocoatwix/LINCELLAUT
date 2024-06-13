import time

startTime = time.time()

Z_N = IntegerModRing(5)
ZP_N, lamb = Z_N["λ"].objgen()

goodPoly = 3 + 3*lamb + 4*lamb^2 + 4*lamb^3 + 1*lamb^4
tempPoly = 0
order = 1

#Find the order
while (True):
    tempPoly = (tempPoly + 1)*lamb - 1
    if (tempPoly % goodPoly == 0):
        print("Order:", order)
        break
        
    order += 1

endTime = time.time()
print("Runtime:", endTime - startTime, "seconds")
