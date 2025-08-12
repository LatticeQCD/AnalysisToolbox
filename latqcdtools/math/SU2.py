# 
# SU2.py                                                               
#
# D. Clarke
#

import numpy as np
from latqcdtools.math.math import SUN, log
from latqcdtools.base.check import checkType, checkDomain


# Pauli matrices
sigma = {
    1: np.array([[0, 1],
                 [1, 0]]),
    2: np.array([[0, -1j],
                 [1j, 0]]),
    3: np.array([[1, 0],
                 [0, -1]])
}


class SU2(SUN):

    def __new__(cls, mat=None):
        return super(SU2, cls).__new__(cls, N=2, mat=mat)

    def __repr__(self) -> str:
        return "SU(2)"

    def projectPauli(self,a) -> float:
        checkType('int',a=a)
        checkDomain(a,[1,2,3])
        M = -1j*log(self)
        return 0.5 * np.trace(sigma[a] @ M).real

    def setToRandom(self):
        u0, u1, u2, u3 = np.random.randn(4)
        norm = np.sqrt(u0**2 + u1**2 + u2**2 + u3**2)
        u0, u1, u2, u3 = u0/norm, u1/norm, u2/norm, u3/norm
        Urand = np.array([[ u0 + 1j * u3, u1 + 1j * u2],
                          [-u1 + 1j * u2, u0 - 1j * u3]])
        self.setToMatrix(Urand)
