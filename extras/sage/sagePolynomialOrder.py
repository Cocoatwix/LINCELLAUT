import time

startTime = time.time()

Z_N = IntegerModRing(13)
ZP_N, lamb = Z_N["λ"].objgen()

goodPoly = 8 + 1*lamb + 11*lamb^2 + 4*lamb^3 + 12*lamb^5
order = 1

#Find the order
while (True):
    tempPoly = lamb^order - 1
    if (tempPoly % goodPoly == 0):
        print("Order:", order)
        break
        
    order += 1

endTime = time.time()
print("Runtime:", endTime - startTime, "seconds")
