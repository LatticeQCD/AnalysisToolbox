i#
# diracFreespectra.py
#
# J. Goswami
#
# Python implementation of the spectrum of the free Lattice QCD dirac operators.
#

import numpy as np
import itertools
import latqcdtools.base.logger as logger
from latqcdtools.math.math import id

logger.warn("This might be broken... Use with caution.")


def _momenta_1d(L: int, antiperiodic: bool = False) -> np.ndarray:
    """
    Return exactly L momentum modes.
    n = -L//2, ..., L//2-1
    p = 2π (n + shift)/L, shift=1/2 for antiperiodic
    Works for even/odd L.
    """
    n = np.arange(-L // 2, L // 2, dtype=int)
    shift = 0.5 if antiperiodic else 0.0
    return 2.0 * np.pi * (n + shift) / L


class GammaMatrix:
    """The 4x4 gamma matrices used in Euclidean QFT."""

    def __repr__(self) -> str:
        return "GammaMatrix()"

    def g(self, i=1):
        if i == 1:
            gamma = np.array([[0.0, 0.0, 0.0, -1.0j],
                              [0.0, 0.0, -1.0j, 0.0],
                              [0.0, 1.0j, 0.0, 0.0],
                              [1.0j, 0.0, 0.0, 0.0]], dtype=complex)
        elif i == 2:
            gamma = np.array([[0.0, 0.0, 0.0, -1.0],
                              [0.0, 0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0, 0.0],
                              [-1.0, 0.0, 0.0, 0.0]], dtype=complex)
        elif i == 3:
            gamma = np.array([[0.0, 0.0, -1.0j, 0.0],
                              [0.0, 0.0, 0.0, 1.0j],
                              [1.0j, 0.0, 0.0, 0.0],
                              [0.0, -1.0j, 0.0, 0.0]], dtype=complex)
        elif i == 4:
            gamma = np.array([[0.0, 0.0, -1.0, 0.0],
                              [0.0, 0.0, 0.0, -1.0],
                              [-1.0, 0.0, 0.0, 0.0],
                              [0.0, -1.0, 0.0, 0.0]], dtype=complex)
        else:
            raise ValueError("Gamma index must be 1..4")
        return gamma

    def g5(self):
        gamma = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0],
                          [0.0, 0.0, -1.0, 0.0],
                          [0.0, 0.0, 0.0, -1.0]], dtype=complex)
        return gamma


class DiracOp(GammaMatrix):

    def __repr__(self) -> str:
        return f"DiracOp(L=({self.Lx},{self.Ly},{self.Lz},{self.Lt}), fermion='{self.fermion}')"

    def __init__(self, Lx, Ly, Lz, Lt, fermion="Wilson", bc_t="anti"):
        """
        Represent Dirac operator on a Euclidean spacetime lattice.

        bc_t: 'anti' (default) or 'peri'
        """
        self.Lx, self.Ly, self.Lz, self.Lt = Lx, Ly, Lz, Lt
        self.fermion = fermion
        if bc_t not in ("anti", "peri"):
            raise ValueError("bc_t must be 'anti' or 'peri'")
        self.bc_t = bc_t

        # cache these (avoids rebuilding each call)
        self._gammas = np.array([self.g(1), self.g(2), self.g(3), self.g(4)], dtype=complex)
        self._I4 = id(4).astype(complex)

    def p(self):
        """Momentum values px,py,pz,pt. Spatial periodic; temporal controlled by bc_t."""
        px = _momenta_1d(self.Lx, antiperiodic=False)
        py = _momenta_1d(self.Ly, antiperiodic=False)
        pz = _momenta_1d(self.Lz, antiperiodic=False)
        pt = _momenta_1d(self.Lt, antiperiodic=(self.bc_t == "anti"))
        return px, py, pz, pt

    def WilsonOp(self, p, mass):
        """Free Wilson operator with r=1."""
        p = np.asarray(p, dtype=float)
        if p.shape != (4,):
            raise ValueError("p must be length 4")

        D = mass * self._I4
        for mu in range(4):
            D = D + 1j * np.sin(p[mu]) * self._gammas[mu] + (1.0 - np.cos(p[mu])) * self._I4
        return D

    def DWMobius4D(self, p, mass, M=1, b=1.5, c=0.5, Ls=12):
        """4D effective Möbius domain-wall operator."""
        p = np.asarray(p, dtype=float)
        if p.shape != (4,):
            raise ValueError("p must be length 4")

        gamma_5 = self.g5()
        I4 = self._I4

        hkNum = (b + c) * self.WilsonOp(p, -M)
        hkDen = 2 * I4 + (b - c) * self.WilsonOp(p, -M)

        # Avoid explicit inverse:
        hk = gamma_5 @ np.linalg.solve(hkDen, hkNum)

        A = np.linalg.matrix_power(I4 + hk, Ls)
        B = np.linalg.matrix_power(I4 - hk, Ls)

        sgnHk = np.linalg.solve(A + B, A - B)

        Ddw = (1 + mass) / 2 * I4 + (1 - mass) / 2 * (gamma_5 @ sgnHk)
        return Ddw

    def eigvalues(self, mass, M=1, b=1.5, c=0.5, Ls=12, flatten=False):
        """Return eigenvalues over all lattice momenta. Shape: (Nmom,4) or flat (4*Nmom,)"""
        eig = []
        px, py, pz, pt = self.p()

        for p4 in itertools.product(px, py, pz, pt):
            p4 = np.array(p4, dtype=float)
            if self.fermion == "Wilson":
                D = self.WilsonOp(p4, mass)
            elif self.fermion == "DwMobius":
                D = self.DWMobius4D(p4, mass, M, b, c, Ls)
            else:
                logger.warn("Unknown fermion type. Falling back to Wilson.")
                D = self.WilsonOp(p4, mass)

            eig.append(np.linalg.eigvals(D))

        eig = np.array(eig, dtype=complex)
        return eig.ravel() if flatten else eig
