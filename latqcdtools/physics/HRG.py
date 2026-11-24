#
# HRG.py
#
# J. Goswami, D. Clarke, Sipaz Sharma
#
# Collection of methods pertaining to hadron resonance gas calculations. Based in part 
# on code by F. Karsch.
#

import numpy as np
import sympy as sy
import warnings
import functools
from mpmath import mp
from scipy.special import lambertw, kn
from sympy import Sum, symbols, Indexed, lambdify, LambertW, exp
import latqcdtools.base.logger as logger
from latqcdtools.math.num_int import integrateFunction
from latqcdtools.base.check import ignoreUnderflow, checkType
from latqcdtools.base.utilities import envector, unvector


ignoreUnderflow()
warnings.filterwarnings("error")
warnings.filterwarnings("ignore", category=DeprecationWarning)


APPROX_KAON_MASS = 500   # Mass cutoff in [MeV] for Boltzmann approximation.


def RMS_mass(Nt, T):
    a = 6924.46
    b = 31881.4
    c = -1.02357e+06
    mkaon = 493.68
    mpion = 140.0
    x = 1.0/Nt/T*197.3
    pion_rms_mass = mpion + a*x**2+b*x**4+c*x**6
    del2 = pion_rms_mass**2 - mpion**2
    kaon_rms_mass = np.sqrt(mkaon**2+del2)
    return pion_rms_mass, kaon_rms_mass


def LCP_init_NS0(muB):
    """ 
    Give a good initial guess for NS=0 LCP. 
    """
    dS = 0.214
    eS = 0.161
    dQ = 0.0211
    eQ = 0.106
    muQ = -dQ / (1.0 + eQ * muB)
    muS = dS / (1.0 + eS * muB)
    return muQ, muS


def dmuh(order, muh):
    """ 
    d^order/dmuh^order derivative of muh 
    """
    if order == 0:
        return muh
    elif order == 1:
        return 1
    else:
        return 0


def _checkOrders(B,Q,S,C,I):
    """
    Want to avoid simultaneous muQ and muI derivatives.
    """
    checkType("int",B=B)
    checkType("int",Q=Q)
    checkType("int",S=S)
    checkType("int",C=C)
    checkType("int",I=I)
    if (Q>0) and (I>0):
        logger.TBRaise('Do not allow muQ and muI derivatives simulataneously')


def _checkMus(muB,muQ,muS,muC,muI):
    """
    Want to avoid simultaneous muQ and muI. 
    """
    if (muQ!=0) and (muI!=0):
        logger.TBRaise('Do not allow muQ and muI simulataneously')
    

class HRGbase:

    """ 
    Hadron resonance gas base class. Here we collect methods and attributes that all HRG-type classes should have
    in common. Mass=mass of the Hadron/resonance , g=spin degenerecy , w= fermi(-1)/bose(1) statistics.
    B, Q, S, and C are respectively the baryon number, electric charge, strangeness, and charm of each state. 
    """

    def __init__(self, Mass, g, w, B, S, Q, C=None):
        self.Mass = Mass
        self.g  = g
        self.w  = w
        self.B  = B
        self.Q  = Q
        self.S  = S
        self.I3 = Q-(B+S)/2 # Gell-Mann-Nishijima 
        if C is None:
            self.C = np.zeros(len(Mass))
        else:
            self.C = C

    def __repr__(self) -> str:
        return "HRGbase"

    def muN_div_T(self, k, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, muI_div_T):
        """ 
        mu_X * N_X, X = (B,Q,S,C,I) 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        return self.B[k]*muB_div_T + self.Q[k]*muQ_div_T + self.S[k]*muS_div_T + self.C[k]*muC_div_T

    def z(self, k, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, muI_div_T):
        """ 
        e^(mu_X*N_X/T) , X = (B,Q,S,C,I) 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        return np.exp(self.muN_div_T(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T))


class HRG(HRGbase):
    """ 
    HRG implemented through Taylor expasion of logarithm. For more information please see, e.g.
    Physics Letters B 695 (2011) 136–142 or especially arXiv:2011.02812.
    You can optionally adjust NMAX_light and NMAX_heavy, which control the number of terms to keep in the
    Taylor expansion for species that are respectively lighter and heavier than the Kaon.
    Our pressure is given in terms of a Taylor series involving modified Bessel functions of the second kind, which
    needs to be truncated at some order. These functions get strongly suppressed when their argument is large.
    In our case, this argument is proportional to the mass. Hence we will sometimes take the Boltzmann
    approximation, i.e. that the mass is large compared to the temperature. In this limit, even fewer terms of the
    expansion need to be kept. Doing so boosts performance.

    We work in the grand canonical ensemble, i.e. we consider P(V,T,mu_i/T). Hence in this class, derivatives w.r.t.
    one of those N_chemical_potentials + 2 variables assume all others are held fixed. 
    """

    def __init__(self, Mass, g, w, B, S, Q, C=None, NMAX_light=21, NMAX_heavy=2):
        HRGbase.__init__(self, Mass, g, w, B, S, Q, C)
        self.NMAX_light = NMAX_light
        self.NMAX_heavy = NMAX_heavy

    def __repr__(self) -> str:
        return "HRG"

    # Numerical backend. HRG works in double precision; HRG_mp overrides these to work in arbitrary precision.
    _PI = np.pi

    def _asarray(self, a):
        """ 
        Convert input to an array of the backend number type. 
        """
        return np.asarray(a, dtype=float)

    def _zeros(self, shape):
        """ 
        Array of zeros of the backend number type. 
        """
        return np.zeros(shape)

    def _exp(self, x):
        """ 
        Elementwise exponential of the backend. 
        """
        return np.exp(x)

    def _kn(self, nu, x):
        """ 
        Elementwise modified Bessel function of the second kind K_nu(x), integer nu, of the backend. 
        """
        return kn(nu, x)

    def Nmax(self, k):
        """ Keep Nmax terms of the Taylor expansion. """
        if self.Mass[k] < APPROX_KAON_MASS:
            return self.NMAX_light
        else:
            return self.NMAX_heavy

    # n represents the nth order of a Taylor expansion of a logarithm.
    # k represents the kth state that appears in the table of resonances.
    # n=1 is Boltzmann approximation

    def factor(self, k, n, T):
        """ 
        (m/T)^2 g eta^(n+1) / 2pi^2 n^2 
        """
        return (self.Mass[k]/T)**2 * self.g[k] * self.w[k]**(n+1) / (2*np.pi**2*n**2)

    def _stateSum(self, T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel):
        """ 
        Sum over all states k and all kept orders n of the logarithm expansion of
            weight(n) * factor(k,n,T) * z_k^n * kernel(n, m_k, x),   x = n m_k/T,
        vectorized over states, temperatures, and chemical potentials. The sum over n runs up to Nmax-1, as before.

        Args:
            T (float or array-like): temperature.
            muB_div_T (float or array-like): muB/T.
            muQ_div_T (float or array-like): muQ/T.
            muS_div_T (float or array-like): muS/T.
            muC_div_T (float or array-like): muC/T.
            weight (func): charge weights of all states at order n, called as weight(n), returning an array
                over states.
            kernel (func): called as kernel(n, m, x, T), with m and x broadcast over states and the input shape,
                and T broadcast to the input shape.

        Returns:
            float or array-like: the sum, with the broadcast shape of T and the chemical potentials.
        """
        T, muB, muQ, muS, muC = np.broadcast_arrays(*[self._asarray(a) for a in
                                                      (T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T)])
        shape  = T.shape
        expand = (slice(None),) + (np.newaxis,)*len(shape)
        light  = np.asarray(self.Mass, dtype=float) < APPROX_KAON_MASS
        Mass   = self._asarray(self.Mass)
        g      = np.asarray(self.g, dtype=float)
        w      = np.asarray(self.w, dtype=float)
        result = self._zeros(shape)
        for group, nmax in ((light, self.NMAX_light), (~light, self.NMAX_heavy)):
            if not np.any(group):
                continue
            m      = Mass[group][expand]
            expo   = ( np.asarray(self.B, dtype=float)[group][expand]*muB + np.asarray(self.Q, dtype=float)[group][expand]*muQ
                     + np.asarray(self.S, dtype=float)[group][expand]*muS + np.asarray(self.C, dtype=float)[group][expand]*muC )
            for n in range(1, nmax):
                x       = n*m/T
                factor  = (m/T)**2 * g[group][expand] * w[group][expand]**(n+1) / (2*self._PI**2*n**2)
                result += np.sum( weight(n)[group][expand] * factor * self._exp(n*expo) * kernel(n, m, x, T), axis=0 )
        return result[()]

    def _chargeWeight(self, B_order, S_order, Q_order, C_order, I_order):
        """ 
        Charge weights (B n)^B_order (S n)^S_order (Q n)^Q_order (C n)^C_order (I3 n)^I_order as a function of n. 
        """
        B  = np.asarray(self.B, dtype=float)
        S  = np.asarray(self.S, dtype=float)
        Q  = np.asarray(self.Q, dtype=float)
        C  = np.asarray(self.C, dtype=float)
        I3 = np.asarray(self.I3, dtype=float)
        def weight(n):
            return (B*n)**B_order * (S*n)**S_order * (Q*n)**Q_order * (C*n)**C_order * (I3*n)**I_order
        return weight

    def P_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        def weight(n):
            return np.ones(len(self.Mass))
        def kernel(n, m, x, T):
            return self._kn(2, x)
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def E_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        def weight(n):
            return np.ones(len(self.Mass))
        def kernel(n, m, x, T):
            return self._kn(2, x)*3 + self._kn(1, x)*x
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def theta_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        Trace anomaly
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        return self.E_div_T4(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T) - 3*self.P_div_T4(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)

    def S_div_T3(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        s = e + p - mu_i n_i 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        NB = self.gen_chi(T, B_order=1, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NS = self.gen_chi(T, B_order=0, S_order=1, Q_order=0, C_order=0, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NQ = self.gen_chi(T, B_order=0, S_order=0, Q_order=1, C_order=0, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NC = self.gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=1, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NI = self.gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=1, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        muxN_div_T = NB*muB_div_T + NS*muS_div_T + NQ*muQ_div_T + NC*muC_div_T + NI*muI_div_T
        return self.E_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T) \
                + self.P_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T,muQ_div_T=muQ_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T) - muxN_div_T

    def ddT_E_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        d(E/T^4)/dT at fixed mu/T
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        def weight(n):
            return np.ones(len(self.Mass))
        def kernel(n, m, x, T):
            return m*n*(self._kn(0, x)*n*m + self._kn(1, x)*T)/T**3
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def ddT_P_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        d(P/T^4)/dT at fixed mu/T
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        def weight(n):
            return np.ones(len(self.Mass))
        def kernel(n, m, x, T):
            return m*n*self._kn(1, x)/T**2
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def CV_div_T3_mu0(self, T):
        """ 
        C_V/T^3 at mu=0. 
        """
        return 4*self.E_div_T4(T, 0, 0, 0, 0, 0) + T*self.ddT_E_div_T4(T, 0, 0, 0, 0, 0)

    def ddT_S_div_T3(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        d(s/T^3)/dT at fixed mu/T 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        ddTNB = self.ddT_gen_chi(T, B_order=1, S_order=0, Q_order=0, C_order=0, I_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        ddTNS = self.ddT_gen_chi(T, B_order=0, S_order=1, Q_order=0, C_order=0, I_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        ddTNQ = self.ddT_gen_chi(T, B_order=0, S_order=0, Q_order=1, C_order=0, I_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        ddTNC = self.ddT_gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=1, I_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        ddTNI = self.ddT_gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=1,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        muxddTN_div_T = ddTNB*muB_div_T + ddTNS*muS_div_T + ddTNQ*muQ_div_T + ddTNC*muC_div_T + ddTNI*muI_div_T
        return self.ddT_E_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T) \
                + self.ddT_P_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T) - muxddTN_div_T

    def gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        chi_BQSC
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        weight = self._chargeWeight(B_order, S_order, Q_order, C_order, I_order)
        def kernel(n, m, x, T):
            return self._kn(2, x)
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def ddT_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        d(chi_BQSC)/dT at fixed mu/T
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        weight = self._chargeWeight(B_order, S_order, Q_order, C_order, I_order)
        def kernel(n, m, x, T):
            return m*n*self._kn(1, x)/T**2
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def d2dT2_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        d^2(chi_BQSC)/dT^2 at fixed mu/T
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        weight = self._chargeWeight(B_order, S_order, Q_order, C_order, I_order)
        def kernel(n, m, x, T):
            return m*n*(m*n*self._kn(0, x)/T**4 - 3*self._kn(1, x)/T**3)
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def gen_ddmuh_E_div_T4(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        Arbitrary mu/T derivatives of E/T^4 at fixed T
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        weight = self._chargeWeight(B_order, S_order, Q_order, C_order, I_order)
        def kernel(n, m, x, T):
            return self._kn(2, x)*3 + self._kn(1, x)*x
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def gen_ddmuh_P_div_T4(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        Arbitrary mu/T derivatives of P/T^4 at fixed T
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        weight = self._chargeWeight(B_order, S_order, Q_order, C_order, I_order)
        def kernel(n, m, x, T):
            return self._kn(2, x)
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)

    def gen_ddmuh_S_div_T3(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        s = e + p - mu_i n_i 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        NB = self.gen_chi(T, B_order=1, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NS = self.gen_chi(T, B_order=0, S_order=1, Q_order=0, C_order=0, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NQ = self.gen_chi(T, B_order=0, S_order=0, Q_order=1, C_order=0, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NC = self.gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=1, I_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        NI = self.gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=1, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        dmuhNB = self.gen_chi(T, B_order=1+B_order, S_order=S_order, Q_order=Q_order, C_order=C_order, I_order=I_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        dmuhNS = self.gen_chi(T, B_order=B_order, S_order=1+S_order, Q_order=Q_order, C_order=C_order, I_order=I_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        dmuhNQ = self.gen_chi(T, B_order=B_order, S_order=S_order, Q_order=1+Q_order, C_order=C_order, I_order=I_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        dmuhNC = self.gen_chi(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=1+C_order, I_order=I_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        dmuhNI = self.gen_chi(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=C_order, I_order=1+I_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)
        return self.gen_ddmuh_E_div_T4(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=C_order, I_order=I_order, 
                                       muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T) \
                + self.gen_ddmuh_P_div_T4(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=C_order, I_order=I_order,
                                          muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T) \
                - muB_div_T*dmuhNB - muQ_div_T*dmuhNQ - muS_div_T*dmuhNS - muC_div_T*dmuhNC - muI_div_T*dmuhNI \
                - NB*dmuh(B_order, muB_div_T) - NQ*dmuh(Q_order, muQ_div_T) - NS*dmuh(S_order, muS_div_T) - NC*dmuh(C_order, muC_div_T) - NI*dmuh(I_order,muI_div_T)

    def gen_chi_RMS(self, T, Nt, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0., muI_div_T=0.):
        # rms_mass[0] is for pions and rms_mass[1] is for kaons
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        _checkOrders(B=B_order,Q=Q_order,S=S_order,C=C_order,I=I_order)
        rms_mass = RMS_mass(Nt, T)
        chi = 0.0
        for k in range(len(self.Mass)):
            if 140 >= self.Mass[k] >= 130:
                for n in range(1, 20):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n )**S_order \
                                                  * (self.Q[k]*n )**Q_order \
                                                  * (self.C[k]*n )**C_order \
                                                  * (self.I3[k]*n)**I_order \
                                                  * self.w[k]**(n+1) * self.g[k] * (rms_mass[0]/T)**2 \
                                                  * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)**n \
                                                  * self._kn(2, (n*rms_mass[0]/T)) / (np.pi*n)**2 / 2
            elif 500 >= self.Mass[k] >= 490:
                for n in range(1, 10):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n )**S_order \
                                                  * (self.Q[k]*n )**Q_order \
                                                  * (self.C[k]*n )**C_order \
                                                  * (self.I3[k]*n)**I_order \
                                                  * self.w[k]**(n+1) * self.g[k] * (rms_mass[1]/T)**2 \
                                                  * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)**n \
                                                  * self._kn(2, (n*rms_mass[1]/T)) / (np.pi*n)**2 / 2
            else:
                for n in range(1, 2):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n )**S_order \
                                                  * (self.Q[k]*n )**Q_order \
                                                  * (self.C[k]*n )**C_order \
                                                  * (self.I3[k]*n)**I_order \
                                                  * self.factor(k, n, T) \
                                                  * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T)**n \
                                                  * self._kn(2, (n*self.Mass[k]/T))
        return chi

    def genChiFlavor(self, T, u_order=0, d_order=0, s_order=0, c_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0., muI_div_T=0.):
        """ 
        chi_udsc 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        logger.warn('HRGFlavor basis is not yet tested for finite mu...')
        B = np.asarray(self.B, dtype=float)
        Q = np.asarray(self.Q, dtype=float)
        S = np.asarray(self.S, dtype=float)
        C = np.asarray(self.C, dtype=float)
        def weight(n):
            return ((B+Q-C)*n)**u_order * ((2.0*B-Q+S)*n)**d_order * (-S*n)**s_order * (C*n)**c_order
        def kernel(n, m, x, T):
            return self._kn(2, x)
        return self._stateSum(T, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, weight, kernel)
# TODO: gen_chi can be implemented using multi-dimensional Leibniz rule


DEFAULT_DPS = 50   # Default number of decimal digits for HRG_mp.


def suggestedDps(T, deltaM, digits=16):
    """ 
    Suggest a working precision for HRG_mp. The exponent range of mpmath numbers is unbounded, so Boltzmann
    factors never underflow, whatever the precision. Precision is needed only when a result comes out of a
    cancellation. If terms of size A cancel to leave a result of size A exp(-deltaM/T), about deltaM/(T ln 10)
    digits are lost. Example: along a strangeness-neutral line with fixed nQ/nB at low T, nQ is a pi+/pi-
    cancellation that must leave something of the size of the nucleon contribution, so deltaM ~ m_N - m_pi. 

    Args:
        T (float): temperature in [MeV].
        deltaM (float): mass scale in [MeV] controlling the cancellation.
        digits (int, optional): digits wanted in the result. Defaults to 16.

    Returns:
        int: suggested number of decimal digits.
    """
    return int(np.ceil(deltaM/(T*np.log(10)))) + digits


@functools.lru_cache(maxsize=2**16)
def _besselK_mp(nu, x, prec):
    """ 
    Cached mpmath K_nu(x) at binary precision prec, for integer nu >= 0. Orders above 1 use the recurrence
    K_nu(x) = K_{nu-2}(x) + 2(nu-1)/x K_{nu-1}(x), which is exact, adds only positive terms, and saves the
    slow integer-order evaluation in mpmath. 
    """
    with mp.workprec(prec):
        if nu < 2:
            return mp.besselk(nu, x)
        return _besselK_mp(nu-2, x, prec) + 2*(nu-1)/x * _besselK_mp(nu-1, x, prec)


def _besselK_mpScalar(nu, x):
    """ 
    mpmath K_nu(x) at the current precision. 
    """
    return _besselK_mp(int(nu), mp.mpf(x), mp.prec)


def _exp_mpScalar(x):
    """ 
    mpmath exp(x) at the current precision. 
    """
    return mp.exp(x)


def _mpf_scalar(x):
    """ 
    Convert to an mpmath float at the current precision. 
    """
    return mp.mpf(x)


_besselK_mpArray = np.frompyfunc(_besselK_mpScalar, 2, 1)
_exp_mpArray     = np.frompyfunc(_exp_mpScalar, 1, 1)
_mpf_array       = np.frompyfunc(_mpf_scalar, 1, 1)


def _withPrecision(method):
    """ 
    Run an HRG_mp method at the instance's precision. 
    """
    @functools.wraps(method)
    def wrapped(self, *args, **kwargs):
        with mp.workdps(self.dps):
            return method(self, *args, **kwargs)
    return wrapped


class HRG_mp(HRG):
    """ 
    HRG in arbitrary precision using mpmath. Same physics and methods as HRG, but every number is an mpmath
    float, so Boltzmann factors like exp(-m_N/T) never underflow, and cancellations can be resolved by raising
    the precision. Use it when double precision fails, e.g. for T below about 10 MeV. It is much slower than
    HRG. Results are mpmath floats (scalar input) or numpy object arrays of mpmath floats (array input). 

    Each method runs at dps decimal digits, without changing mpmath's global precision. Arithmetic you do with
    the results afterwards happens at the global precision mp.dps, so if your own combinations of results
    involve cancellations, do them inside "with mp.workdps(hrg.dps):". See suggestedDps for choosing dps.

    K_nu(x) is cached across calls, which pays off when many quantities are computed at the same T.
    """

    def __init__(self, Mass, g, w, B, S, Q, C=None, NMAX_light=21, NMAX_heavy=2, dps=DEFAULT_DPS):
        HRG.__init__(self, Mass, g, w, B, S, Q, C, NMAX_light, NMAX_heavy)
        checkType("int", dps=dps)
        self.dps = dps

    def __repr__(self) -> str:
        return "HRG_mp"

    _PI = mp.pi

    def _asarray(self, a):
        return np.array(_mpf_array(np.asarray(a, dtype=object)), dtype=object)

    def _zeros(self, shape):
        return np.full(shape, mp.mpf(0), dtype=object)

    def _exp(self, x):
        return _exp_mpArray(x)

    def _kn(self, nu, x):
        return _besselK_mpArray(nu, x)

    def factor(self, k, n, T):
        """ 
        (m/T)^2 g eta^(n+1) / 2pi^2 n^2 
        """
        return (mp.mpf(self.Mass[k])/T)**2 * self.g[k] * self.w[k]**(n+1) / (2*mp.pi**2*n**2)

    def z(self, k, muB_div_T, muQ_div_T, muS_div_T, muC_div_T, muI_div_T):
        """ 
        e^(mu_X*N_X/T) , X = (B,Q,S,C,I) 
        """
        _checkMus(muB=muB_div_T,muQ=muQ_div_T,muS=muS_div_T,muC=muC_div_T,muI=muI_div_T)
        return mp.exp(self.muN_div_T(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T, muI_div_T=muI_div_T))

    def gen_chi_RMS(self, T, Nt, B_order=0, S_order=0, Q_order=0, C_order=0, I_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0., muI_div_T=0.):
        logger.TBRaise('gen_chi_RMS is not implemented for HRG_mp. Taste-split masses matter near Tpc, where HRG suffices.')


# Every public method of HRG_mp runs at the instance's precision.
for _name in dir(HRG_mp):
    if not _name.startswith('_') and callable(getattr(HRG_mp, _name)):
        setattr(HRG_mp, _name, _withPrecision(getattr(HRG_mp, _name)))



def solveLCP(hrg, T, muB_div_T, rQ, start=None, tol=1e-12, maxiter=100):
    """ 
    Solve NS = 0 and NQ/NB = rQ for muQ/T and muS/T at fixed T and muB/T. Works for HRG and HRG_mp.

    The equations are solved as NS/NB = 0 and NQ/NB - rQ = 0, which are of order one at every T. (Densities
    carry a factor ~exp(-m_N/T), so at low T an absolute criterion on NS and NQ - rQ NB would accept any guess.)
    The Jacobian is exact, since d NX/d(muY/T) are second-order cumulants, and each Newton step solves the 2x2
    system by Cramer's rule, halving the step until the residual decreases. At low T the solutions are
    exponentially small, and NQ comes from a pi+/pi- cancellation, so there you need HRG_mp with enough digits;
    see suggestedDps. Tracking a line down in T works best by using each solution as start for the next T.

    Args:
        hrg (HRG): HRG or HRG_mp object.
        T (float): temperature in [MeV].
        muB_div_T (float): muB/T.
        rQ (float): NQ/NB.
        start (tuple, optional): starting guess (muQ/T, muS/T). Defaults to LCP_init_NS0(muB_div_T).
        tol (float, optional): tolerance on both residuals. Defaults to 1e-12.
        maxiter (int, optional): maximum number of Newton steps. Defaults to 100.

    Raises:
        Exception: no convergence within maxiter steps.

    Returns:
        tuple: muQ/T, muS/T
    """
    checkType("real", rQ=rQ)
    checkType("int", maxiter=maxiter)
    if start is None:
        start = LCP_init_NS0(muB_div_T)
    if isinstance(hrg, HRG_mp):
        with mp.workdps(hrg.dps):
            return _solveLCP(hrg, mp.mpf(T), mp.mpf(muB_div_T), mp.mpf(rQ), mp.mpf(start[0]), mp.mpf(start[1]), tol, maxiter)
    return _solveLCP(hrg, T, muB_div_T, rQ, start[0], start[1], tol, maxiter)


def _LCPresiduals(hrg, T, muB_div_T, muQ_div_T, muS_div_T, rQ):
    """ 
    (NS/NB, NQ/NB - rQ) and the densities NB, NQ, NS. 
    """
    NB = hrg.gen_chi(T, B_order=1, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T)
    NQ = hrg.gen_chi(T, Q_order=1, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T)
    NS = hrg.gen_chi(T, S_order=1, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T)
    return NS/NB, NQ/NB - rQ, NB, NQ, NS


def _solveLCP(hrg, T, muB_div_T, rQ, muQ_div_T, muS_div_T, tol, maxiter):
    """ 
    Newton iteration for solveLCP. 
    """
    resS, resQ, NB, NQ, NS = _LCPresiduals(hrg, T, muB_div_T, muQ_div_T, muS_div_T, rQ)
    for _ in range(maxiter):
        size = max(abs(resS), abs(resQ))
        if size < tol:
            return muQ_div_T, muS_div_T
        mus  = {"muB_div_T": muB_div_T, "muQ_div_T": muQ_div_T, "muS_div_T": muS_div_T}
        X11BQ = hrg.gen_chi(T, B_order=1, Q_order=1, **mus)
        X11BS = hrg.gen_chi(T, B_order=1, S_order=1, **mus)
        X2Q   = hrg.gen_chi(T, Q_order=2, **mus)
        X11QS = hrg.gen_chi(T, Q_order=1, S_order=1, **mus)
        X2S   = hrg.gen_chi(T, S_order=2, **mus)
        # Jacobian of (NS/NB, NQ/NB) w.r.t. (muQ/T, muS/T), using d NX/d(muY/T) = chi_11^XY
        JSQ = (X11QS*NB - NS*X11BQ)/NB**2
        JSS = (X2S  *NB - NS*X11BS)/NB**2
        JQQ = (X2Q  *NB - NQ*X11BQ)/NB**2
        JQS = (X11QS*NB - NQ*X11BS)/NB**2
        det = JSQ*JQS - JSS*JQQ
        dQ  = (-resS*JQS + resQ*JSS)/det
        dS  = (-resQ*JSQ + resS*JQQ)/det
        step = 1
        for _ in range(30):
            newResS, newResQ, newNB, newNQ, newNS = _LCPresiduals(hrg, T, muB_div_T, muQ_div_T + step*dQ, muS_div_T + step*dS, rQ)
            if max(abs(newResS), abs(newResQ)) < size:
                break
            step = step/2
        muQ_div_T, muS_div_T = muQ_div_T + step*dQ, muS_div_T + step*dS
        resS, resQ, NB, NQ, NS = newResS, newResQ, newNB, newNQ, newNS
    if max(abs(resS), abs(resQ)) < tol:
        return muQ_div_T, muS_div_T
    logger.TBRaise(f"LCP not solved at T = {float(T)}, muB/T = {float(muB_div_T)}: NS/NB = {float(resS):.3e}, NQ/NB - rQ = {float(resQ):.3e}")


class EVHRG(HRGbase):

    """ 
    Excluded volume hadron resonance gas. b is excluded volume parameter. 
    """

    def __init__(self, Mass, g, w, B, S, Q, C=None):
        logger.warn("Only implemented for (B,Q,S) basis")
        HRGbase.__init__(self, Mass, g, w, B, S, Q, C)

    def __repr__(self) -> str:
        return "EVHRG"

    def Pid(self, m, g, T):
        return g*(m/T)**2 * kn(2, (m/T)) / np.pi**2 / 2

    def baryon_pressure(self, T, b, Bi, muB_div_T=0.0, muQ_div_T=0.0, muS_div_T=0.0):

        baryon_mass = self.Mass[np.where(self.B == Bi)]
        g_baryon = self.g[np.where(self.B == Bi)]
        X_baryon = self.B[np.where(self.B == Bi)]
        X_charge = self.Q[np.where(self.B == Bi)]
        X_strange = self.S[np.where(self.B == Bi)]

        # Bi=1 Baryon pressure , Bi=-1 for anti baryon pressure
        P = []
        for k in range(len(baryon_mass)):
            P.append(self.Pid(baryon_mass[k], g_baryon[k], T))

        P = np.array(P)

        mB, mQ, mS, i, ci, bi, qi, si, temp = symbols(
            'mB, mQ , mS, i, ci, bi, qi, si, temp')

        F_pressure = (1 / (b*(T/197.3)**3)) * LambertW(
            (b*(T/197.3)**3) * Sum(Indexed('ci', i) * exp(mB * Indexed('bi', i) + mQ * Indexed('qi', i) + mS * Indexed('si', i)),
                                   (i, 0, len(P) - 1)))
        f = lambdify((mB, mQ, mS, ci, bi, qi, si, temp), F_pressure,
                     modules=['scipy', {'LambertW': lambertw}])

        pressure = f(muB_div_T/T, muQ_div_T/T, muS_div_T/T,
                     P, X_baryon, X_charge, X_strange, T).real

        return pressure

    def gen_chi(self, T, b, Bi, B_order=0, Q_order=0, S_order=0, muB_div_T=0.0, muQ_div_T=0.0, muS_div_T=0.0):

        baryon_mass = self.Mass[np.where(self.B == Bi)]
        g_baryon = self.g[np.where(self.B == Bi)]
        X_baryon = self.B[np.where(self.B == Bi)]
        X_charge = self.Q[np.where(self.B == Bi)]
        X_strange = self.S[np.where(self.B == Bi)]

        # Bi=1 Baryon pressure , Bi=-1 for anti baryon pressure
        P = []
        for k in range(len(baryon_mass)):
            P.append(self.Pid(baryon_mass[k], g_baryon[k], T))

        P = np.array(P)

        # ci = ideal gas pressure of individual particles
        mB, mQ, mS, i, ci, bi, qi, si, temp, be = symbols(
            'mB, mQ , mS, i, ci, bi, qi, si, temp, be')

        F_pressure = (1 / be) * LambertW(
            be * Sum(
                Indexed('ci', i) * exp(mB * Indexed('bi', i) + mQ *
                                       Indexed('qi', i) + mS * Indexed('si', i)),
                (i, 0, len(P) - 1)))
        expr_chi = sy.diff(F_pressure, mB, B_order, mQ, Q_order, mS, S_order)

        f = lambdify((mB, mQ, mS, ci, bi, qi, si, temp, be),
                     expr_chi, modules=['scipy', {'LambertW': lambertw}])

        chi_num = f(muB_div_T/T, muQ_div_T/T, muS_div_T/T, P,
                    X_baryon, X_charge, X_strange, T, b*(T/197.3)**3).real

        return chi_num
