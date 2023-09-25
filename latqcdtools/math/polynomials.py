# 
# polynomials.py                                                               
# 
# J. Goswami
# 
# A simple class for constructing polynomial-type and ratioal-type functions. You can construct a polynomial of only
# even powers up to fourth order using, for example,
#   p = Polynomial([A0, 0., A2, 0. A4])
#   p(x)
#

from latqcdtools.base.check import checkType

class Polynomial:

    def __init__(self, coeffs=None):
        checkType(coeffs,"array")
        self.coeffs = coeffs

    @property
    def __repr__(self) -> str:
        return "Polynomial{0}".format((str(self.coeffs)))

    def __call__(self, x):
        res = 0.0
        for index_num in range(len(self.coeffs)):
            res += self.coeffs[index_num] *x**index_num
        return res


class Rational:

    def __init__(self, num_coeffs, den_coeffs):
        checkType(num_coeffs,"array")
        checkType(den_coeffs,"array")
        self.num_coeffs = num_coeffs
        self.den_coeffs = den_coeffs

    @property
    def __repr__(self) -> str:
        return "Rational function" + (str(self.num_coeffs) + str(self.den_coeffs))

    def __call__(self, x):
        res = 0.0
        res_den = 0.0
        for index_den in range(len(self.den_coeffs)):
            res_den += self.den_coeffs[index_den] *x**index_den
        for index_num in range(len(self.num_coeffs)):
            res += self.num_coeffs[index_num] *x**index_num/res_den
        return res
