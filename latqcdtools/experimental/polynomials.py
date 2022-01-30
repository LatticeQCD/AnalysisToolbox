class Polynomial:
    def __init__(self, coeffs=None):
        self.coeffs = coeffs

    @property
    def __repr__(self):
        return "Polynomial{0}".format((str(self.coeffs)))

    def __call__(self, x):
        res = 0.0
        for index_num in range(len(self.coeffs)):
            res += self.coeffs[index_num] * x ** index_num
        return res


class Rational:

    def __init__(self, num_coeffs, den_coeffs):

        self.num_coeffs = num_coeffs
        self.den_coeffs = den_coeffs

    @property
    def __repr__(self):
        return "Pade" + (str(self.num_coeffs) + str(self.den_coeffs))

    def __call__(self, x):
        res = 0.0
        res_den = 0.0
        for index_den in range(len(self.den_coeffs)):
            res_den += self.den_coeffs[index_den] * x ** index_den
        for index_num in range(len(self.num_coeffs)):
            res += self.num_coeffs[index_num] * x ** index_num / res_den
        return res
