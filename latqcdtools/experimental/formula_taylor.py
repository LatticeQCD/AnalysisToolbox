"""
Calculation of the Higher Order Taylor co-efficients
Create by Jishnu Goswami 18/01/2022
"""

import math as mp
import sympy as sy
from sympy import symbols

order = 6

muB, muQ, muS, s1, s3, s5, s7, s9, q1, q3, q5, q7, q9 = symbols('muB,muQ,muS,s1,s3,s5,s7,s9,q1,q3,q5,q7,q9')

def constraint_expansion(expr):
    tmp = expr.subs({muS: s1 * muB + s3 * muB ** 3 + s5 * muB ** 5 + s7 * muB ** 7 + s9 * muB ** 9,
            muQ: q1 * muB + q3 * muB ** 3 + q5 * muB ** 5 + q7 * muB ** 7 + q9 * muB ** 9})
    return tmp.expand()


P = sum(1 / sy.factorial(i) * 1 / sy.factorial(j) * 1 / sy.factorial(k) * muB ** i * muQ ** j * muS ** k * symbols(
    'X.{}.{}.{}'.format(i, j, k)) for i in range(order) for j in range(order) for k in range(order) if
        (i + j + k) % 2 == 0 and i + j + k < order)

"""
Co-efficients of Pressure
"""
print("Calculations of Co-efficients of Pressure....")
m = 2
F = constraint_expansion(P)
PBm = mp.factorial(m) * F.coeff(muB ** m)

print (PBm)

print("done")

"""
Co-efficients of Number density
"""
print("Calculations of Number density....")

k = m - 1
t = sy.diff(P, muB, 1)

Nb = constraint_expansion(t)

NBk = mp.factorial(k) * Nb.coeff(muB ** k)

print("done")

"""
Co-efficients of chi2B
"""
print("Calculations of chi2B....")

l = m - 2
cexpr = sy.diff(P, muB, 2)

cb = constraint_expansion(cexpr)
CBl = mp.factorial(l) * cb.coeff(muB ** l)


if (l == 0):
    CBl = cb.extract_leading_order(muB)[0][0]
else:
    CBl = mp.factorial(l) * cb.coeff(muB ** l)

print("done")
out_file = open("Taylor_coeffcient_%dorder.txt" % (m), "w")
out_file.write('#Pressure : %d order' % (m) + "\n")
out_file.write(r'P%dB=' % (m))
out_file.write(str(PBm))
out_file.write("\n")
out_file.write('# number density : %d order' % (k) + "\n")
out_file.write(r'NB%dB=' % (k))
out_file.write(str(NBk))
out_file.write("\n")
out_file.write(r'#$\chi^B_2$ : %d order' % (l) + "\n")
out_file.write(r'CB%dB=' % (l))
out_file.write(str(CBl))
out_file.close()
