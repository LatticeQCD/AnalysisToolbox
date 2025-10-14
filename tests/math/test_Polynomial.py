# 
# testPolynomial.py
# 
# D. Clarke
# 
# Testing for the Polynomial and Rational classes.
# 

from latqcdtools.testing import print_results, concludeTest
from latqcdtools.math.polynomials import Polynomial, Rational


def testPolynomial():

    lpass = True

    # 1 + 2 + 4 + 8 = 15
    p = Polynomial([1,1,1,1])
    lpass *= print_results(p(2),15,text="Polynomial test",prec=1e-10)

    r = Rational([1,1,1,1],[1,1,1,1])
    lpass *= print_results(r(2),1,text="Rational test",prec=1e-10)

    concludeTest(lpass)


if __name__ == '__main__':
    testPolynomial()