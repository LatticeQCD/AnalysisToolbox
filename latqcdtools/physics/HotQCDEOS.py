# 
# HotQCDEOS.py                                                               
# 
# J. Goswami 
# 

import numpy as np
import sympy as sy
import warnings
from latqcdtools.math.polynomials import Rational
from latqcdtools.math.optimize import persistentSolve

warnings.filterwarnings("error")
warnings.filterwarnings("ignore", category=DeprecationWarning)

pid = 95 * np.pi ** 2 / 180
x0 = 0.9761
ct = 3.8706

T0 = 154

T = np.arange(130, 305, 5)

a = np.array([0, -8.7704, 3.9200, 0, 0.3419])
b = np.array([0, -1.2600, 0.8425, 0, -0.0475])

r = 0.4


class EOS:
    def __init__(self, temp):
        self.temp = temp

    def __repr__(self) -> str:
        return "EOS"

    def pressure(self):
        """
        calculation of pressure from the parametrization of 1407.6387
        
        """
        T = self.temp / T0
        i, ai, bi, x = sy.symbols('i, ai, bi, x')
        exprP = 0.5 * (1 + sy.tanh(ct * (x - x0))) * (pid + sum(a[i] / x ** i for i in range(1, 5))) / (
                    1 + sum(b[i] / x ** i for i in range(1, 5)))
        pressure = sy.lambdify((ai, bi, x), exprP, modules=[
            'numpy', {'tanh': np.tanh}])(a, b, T)
        return pressure

    def dpressure(self):
        """
        calculation of derivative of pressure from the parametrization of 1407.6387
        """
        T = self.temp / T0
        i, ai, bi, x = sy.symbols('i, ai, bi, x')
        exprP = 0.5 * (1 + sy.tanh(ct * (x - x0))) * (pid + sum(a[i] / x ** i for i in range(1, 5))) \
                / (1 + sum(b[i] / x ** i for i in range(1, 5)))
        exprDerivative = x * sy.diff(exprP, x, 1)
        dPressure = sy.lambdify((ai, bi, x), exprDerivative, modules=[
            'numpy', {'tanh': np.tanh}])(a, b, T)
        return dPressure

    def ObsEoS(self):
        """
        This routine is to calculate pressure , energy density and entropy density for muB = 0
        @return: pressure , energy density and entropy density at muB = 0
        """
        p0 = self.pressure()
        e0 = 3 * p0 + self.dpressure()
        s0 = p0 + e0
        return p0, e0, s0

    def TDerivatives(self, y, pars):
        x = sy.symbols('x')
        a, b, c, d, e, f, g, h, i = sy.symbols('a b c d e f g h i')

        num = a + b * x + c * x ** 2 + d * x ** 3 + e * x ** 4
        den = 1 + f * x + g * x ** 2 + h * x ** 3 + i * x ** 4

        deriv = sy.lambdify((x, a, b, c, d, e, f, g, h, i), sy.diff(num / den, x))

        a0, a1, a2, a3, a4, b1, b2, b3, b4 = pars
        return deriv(y, a0, a1, a2, a3, a4, b1, b2, b3, b4)

    def createCoeffs(self, y, par):
        return Rational([par[0], par[1], par[2], par[3], par[4]], [1, par[5], par[6], par[7], par[8]])(y)

    def ObsEoSfixedmuB(self, Np, qP, muB = 1.0):
        """
        calculation of equation of state( pressure, energy density and
        entropy density for the functional form given in,
        arXiv : 2212.10016 or 1701.04325 [hep-lat]
        @param Np: list of parameters for number density co-efficients
        @param qP: list of parameters for q's
        @param muB: value of the fixed muB
        @return: pressure , energy density and entropy density at fixed muB
        """
        p0 = self.pressure()
        e0 = 3 * p0 + self.dpressure()
        s0 = p0 + e0

        T = T0 / self.temp

        N1B = self.createCoeffs(T, Np[0])
        N3B = self.createCoeffs(T, Np[1])
        N5B = self.createCoeffs(T, Np[2])

        dN1B = (-T) * self.TDerivatives(T, Np[0])
        dN3B = (-T) * self.TDerivatives(T, Np[1])
        dN5B = (-T) * self.TDerivatives(T, Np[2])

        q1 = self.createCoeffs(T, qP[0])
        q3 = self.createCoeffs(T, qP[1])
        q5 = self.createCoeffs(T, qP[2])

        dq1 = (-T) * self.TDerivatives(T, qP[0])
        dq3 = (-T) * self.TDerivatives(T, qP[1])
        dq5 = (-T) * self.TDerivatives(T, qP[2])

        P2B = (N1B + r * q1 * N1B) / 2
        P4B = (N3B + r * (q1 * N3B + 3 * q3 * N1B)) / 4
        P6B = (N5B + r * (q1 * N5B + 3 * q3 * N3B + 5 * q5 * N1B)) / 6

        dP2B = (dN1B + r * (q1 * dN1B + dq1 * N1B)) / 2
        dP4B = (dN3B + r * (dq1 * N3B + q1 * dN3B + 3 * dq3 * N1B + 3 * q3 * dN1B)) / 4
        dP6B = (dN5B + r * (dq1 * N5B + q1 * dN5B + 3 * dq3 * N3B + 3 * q3 * dN3B + 5 * dq5 * N1B + 5 * q5 * dN1B)) / 6

        E2B = (3 * P2B + dP2B - r * dq1 * N1B)
        E4B = (3 * P4B + dP4B - r * (dq1 * N3B + dq3 * N1B))
        E6B = (3 * P6B + dP6B - r * (dq1 * N5B + dq3 * N3B + dq5 * N1B))

        S2B = (4 * P2B + dP2B - N1B - r * ((q1 + dq1) * N1B))
        S4B = (4 * P4B + dP4B - N3B - r * ((q1 + dq1) * N3B + (q3 + dq3) * N1B))
        S6B = (4 * P6B + dP6B - N5B - r * ((q1 + dq1) * N5B + (q3 + dq3) * N3B + (q5 + dq5) * N1B))

        P = p0 + P2B * muB ** 2 + P4B * muB ** 4 + P6B * muB ** 6
        
        NB = N1B * muB + N3B * muB ** 3 + N5B * muB ** 5

        E = e0 + E2B * muB ** 2 + E4B * muB ** 4 + E6B * muB ** 6

        S = s0 + S2B * muB ** 2 + S4B * muB ** 4 + S6B * muB ** 6

        return P, NB, E, S

    def ObsEoSfixedsnB(self, Np, qP, snB = 50):
        """
        calculation of equation of state( pressure,
        energy density and entropy density for the functional
        form given in, arXiv : 2212.10016 or 1701.04325 [hep-lat]
        @param Np: list of parameters for number density co-efficients
        @param qP: list of parameters for q's
        @param snB: For fixed s/nB
        @return: muB/T , pressure , energy density and entropy density at fixed s/nB
        """
        p0 = self.pressure()
        e0 = 3 * p0 + self.dpressure()
        s0 = p0 + e0

        T = T0 / self.temp

        N1B = self.createCoeffs(T, Np[0])
        N3B = self.createCoeffs(T, Np[1])
        N5B = self.createCoeffs(T, Np[2])

        dN1B = (-T) * self.TDerivatives(T, Np[0])
        dN3B = (-T) * self.TDerivatives(T, Np[1])
        dN5B = (-T) * self.TDerivatives(T, Np[2])

        q1 = self.createCoeffs(T, qP[0])
        q3 = self.createCoeffs(T, qP[1])
        q5 = self.createCoeffs(T, qP[2])

        dq1 = (-T) * self.TDerivatives(T, qP[0])
        dq3 = (-T) * self.TDerivatives(T, qP[1])
        dq5 = (-T) * self.TDerivatives(T, qP[2])

        P2B = (N1B + r * q1 * N1B) / 2
        P4B = (N3B + r * (q1 * N3B + 3 * q3 * N1B)) / 4
        P6B = (N5B + r * (q1 * N5B + 3 * q3 * N3B + 5 * q5 * N1B)) / 6

        dP2B = (dN1B + r * (q1 * dN1B + dq1 * N1B)) / 2
        dP4B = (dN3B + r * (dq1 * N3B + q1 * dN3B + 3 * dq3 * N1B + 3 * q3 * dN1B)) / 4
        dP6B = (dN5B + r * (dq1 * N5B + q1 * dN5B + 3 * dq3 * N3B + 3 * q3 * dN3B + 5 * dq5 * N1B + 5 * q5 * dN1B)) / 6

        E2B = (3 * P2B + dP2B - r * dq1 * N1B)
        E4B = (3 * P4B + dP4B - r * (dq1 * N3B + dq3 * N1B))
        E6B = (3 * P6B + dP6B - r * (dq1 * N5B + dq3 * N3B + dq5 * N1B))

        S2B = (4 * P2B + dP2B - N1B - r * ((q1 + dq1) * N1B))
        S4B = (4 * P4B + dP4B - N3B - r * ((q1 + dq1) * N3B + (q3 + dq3) * N1B))
        S6B = (4 * P6B + dP6B - N5B - r * ((q1 + dq1) * N5B + (q3 + dq3) * N3B + (q5 + dq5) * N1B))

        def Psol(x): return p0 + P2B * x ** 2 + P4B * x ** 4 + P6B * x ** 6

        def Esol(x): return e0 + E2B * x ** 2 + E4B * x ** 4 + E6B * x ** 6

        def Ssol(x): return s0 + S2B * x ** 2 + S4B * x ** 4 + S6B * x ** 6

        def NBsol(x): return N1B * x + N3B * x ** 3 + N5B * x ** 5

        def fun(x): return Ssol(x) - snB * NBsol(x)

        xs = np.ones(N1B.size)

        muBsol = persistentSolve(fun, xs)

        return muBsol, Psol(muBsol), NBsol(muBsol), Esol(muBsol), Ssol(muBsol)
