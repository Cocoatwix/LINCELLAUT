
#https://doc.sagemath.org/html/en/reference/rings/sage/rings/quotient_ring.html

p = 31
Z_p = IntegerModRing(p)
Z_pX, lamb = Z_p["λ"].objgen()

Z_pXD.<x> = QuotientRing(Z_pX, Z_pX.ideal(8 + 2*lamb + 4*lamb^3 + 1*lamb^4))

i = 1
while True:
    if (x^i == 1):
        print(i)
        break
    i += 1
