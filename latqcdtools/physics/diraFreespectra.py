# 
# diraFreespectra.py                                                               
# 
# J. Goswami
# 
# Python implementation of the spectrum of the free Lattice QCD dirac operators.
#

import numpy as np
from numpy import ndarray
import latqcdtools.base.logger as logger
import itertools

I4 = np.eye(4)

class GammaMatrix:
    @property
    def __repr__(self):
        return 'Gamma Matrix'

    def g(self, i=1):
        if i == 1:
            gamma: ndarray = np.array([[0.0, 0.0, 0.0, -1.0j],
                                       [0.0, 0.0, -1.0j, 0.0],
                                       [0.0, 1.0j, 0.0, 0.0],
                                       [1.0j, 0.0, 0.0, 0.0]])
        elif i == 2:
            gamma: ndarray = np.array([[0.0, 0.0, 0.0, -1.0],
                                       [0.0, 0.0, 1.0, 0.0],
                                       [0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 0.0]])
        elif i == 3:
            gamma: ndarray = np.array([[0.0, 0.0, -1.0j, 0.0],
                                       [0.0, 0.0, 0.0, 1.0j],
                                       [1.0j, 0.0, 0.0, 0.0],
                                       [0.0, -1.0j, 0.0, 0.0]])
        elif i == 4:
            gamma: ndarray = np.array([[0.0, 0.0, -1.0, 0.0],
                                       [0.0, 0.0, 0.0, -1.0],
                                       [-1.0, 0.0, 0.0, 0.0],
                                       [0.0, -1.0, 0.0, 0.0]])

        return gamma

    def g5(self):
        gamma: ndarray = np.array([[1.0, 0.0, 0.0, 0.0],
                                   [0.0, 1.0, 0.0, 0.0],
                                   [0.0, 0.0, -1.0, 0.0],
                                   [0.0, 0.0, 0.0, -1.0]])
        return gamma


class DiracOp(GammaMatrix):
    @property
    def __repr__(self):
        return "DiracOp"

    def __init__(self, Lx=None, Ly=None, Lz=None, Lt=None, fermion="Wilson"):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.Lt = Lt
        self.fermion = fermion

    def p(self):
        if self.Lt == self.Lx:
            px = 2 * np.pi * np.arange(-self.Lx / 2 + 1, self.Lx / 2 + 1, 1) / self.Lx
            py = 2 * np.pi * np.arange(-self.Ly / 2 + 1, self.Ly / 2 + 1, 1) / self.Ly
            pz = 2 * np.pi * np.arange(-self.Lz / 2 + 1, self.Lz / 2 + 1, 1) / self.Lz
            pt = 2 * np.pi * np.arange(-self.Lt / 2 + 1, self.Lt / 2 + 1, 1) / self.Lt
        else:
            px = 2 * np.pi * np.arange(-self.Lx / 2 + 1, self.Lx / 2 + 1, 1) / self.Lx
            py = 2 * np.pi * np.arange(-self.Ly / 2 + 1, self.Ly / 2 + 1, 1) / self.Ly
            pz = 2 * np.pi * np.arange(-self.Lz / 2 + 1, self.Lz / 2 + 1, 1) / self.Lz
            pt = 2 * np.pi * np.arange(-self.Lt / 2 + 1 + 1 / 2, self.Lt / 2 + 1 + 1 / 2, 1) / self.Lt
        return px, py, pz, pt

    def WilsonOp(self, p, mass):
        gamma = np.array([self.g(1), self.g(2), self.g(3), self.g(4)])
        term = sum(1j * np.sin(p[i]) * gamma[i] +
                   (1 - np.cos(p[i])) * I4 for i in range(len(gamma)))
        massterm = mass * I4
        DWilson = term + massterm
        return DWilson

    def DWMobius4D(self, p, mass, M=1, b=1.5, c=0.5, Ls=12):
        gamma_5 = self.g5()
        hkNum = (b + c) * self.WilsonOp(p, -M)
        hkDen = 2 * I4 + (b - c) * self.WilsonOp(p, -M)
        hk = gamma_5 @ hkNum @ np.linalg.inv(hkDen)
        sgnHk = np.matmul((np.linalg.matrix_power(I4 + hk, Ls) - np.linalg.matrix_power(I4 - hk, Ls)), \
                          np.linalg.inv(np.linalg.matrix_power(I4 + hk, Ls) + np.linalg.matrix_power(I4 - hk, Ls)))
        Ddw = (1 + mass) / 2 * I4 + (1 - mass) / 2 * np.matmul(gamma_5, sgnHk)
        return Ddw

    def eigvalues(self, mass, M=1, b=1.5, c=0.5, Ls=12):
        eig = []
        px, py, pz, pt = self.p()
        if self.fermion == "Wilson":
            for (pi, pj, pk, pl) in itertools.product(px, py, pz, pt):
                D = self.WilsonOp([pi, pj, pk, pl], mass)
                eig.append(np.linalg.eig(D)[0])
        elif self.fermion == 'DwMobius':
            for (pi, pj, pk, pl) in itertools.product(px, py, pz, pt):
                D = self.DWMobius4D([pi, pj, pk, pl], mass, M, b, c, Ls)
                eig.append(np.linalg.eig(D)[0])
        else:
            logger.warn('Fermion Not implemented yet, here is spectrum for Wilson fermion : option : Wilson , '
                            'DwMobius')
            for (pi, pj, pk, pl) in itertools.product(px, py, pz, pt):
                D = self.WilsonOp([pi, pj, pk, pl], mass)
                eig.append(np.linalg.eig(D)[0])

        return np.array(eig)
