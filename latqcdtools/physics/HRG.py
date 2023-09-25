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
from scipy.special import lambertw, kn
from sympy import Sum, symbols, Indexed, lambdify, LambertW, exp
import latqcdtools.base.logger as logger
from latqcdtools.math.num_int import integrateFunction
from latqcdtools.base.check import ignoreUnderflow
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
    """ Give a good initial guess for NS=0 LCP. """
    dS = 0.214
    eS = 0.161
    dQ = 0.0211
    eQ = 0.106
    muQ = -dQ / (1.0 + eQ * muB)
    muS = dS / (1.0 + eS * muB)
    return muQ, muS


def dmuh(order, muh):
    """ d^order/dmuh^order derivative of muh """
    if order == 0:
        return muh
    elif order == 1:
        return 1
    else:
        return 0


class HRGbase:

    """ Hadron resonance gas base class. Here we collect methods and attributes that all HRG-type classes should have
        in common. Mass=mass of the Hadron/resonance , g=spin degenerecy , w= fermi(-1)/bose(1) statistics.
        B, Q, S, and C are respectively the baryon number, electric charge, strangeness, and charm of each state. """

    def __init__(self, Mass, g, w, B, S, Q, C=None):
        self.Mass = Mass
        self.g = g
        self.w = w
        self.B = B
        self.Q = Q
        self.S = S
        if C is None:
            self.C = np.zeros(len(Mass))
        else:
            self.C = C

    def __repr__(self) -> str:
        return "HRGbase"

    def muN_div_T(self, k, muB_div_T, muQ_div_T, muS_div_T, muC_div_T):
        """ mu_X * N_X, X = (B,Q,S,C) """
        return self.B[k]*muB_div_T + self.Q[k]*muQ_div_T + self.S[k]*muS_div_T + self.C[k]*muC_div_T

    def z(self, k, muB_div_T, muQ_div_T, muS_div_T, muC_div_T):
        """ e^(mu_X*N_X/T) , X = (B,Q,S,C) """
        return np.exp(self.muN_div_T(k, muB_div_T, muQ_div_T, muS_div_T, muC_div_T))


class HRG(HRGbase):
    """ HRG implemented through Taylor expasion of logarithm. For more information please see, e.g.
        Physics Letters B 695 (2011) 136–142 or especially arXiv:2011.02812.
        You can optionally adjust NMAX_light and NMAX_heavy, which control the number of terms to keep in the
        Taylor expansion for species that are respectively lighter and heavier than the Kaon.
        Our pressure is given in terms of a Taylor series involving modified Bessel functions of the second kind, which
        needs to be truncated at some order. These functions get strongly suppressed when their argument is large.
        In our case, this argument is proportional to the mass. Hence we will sometimes take the Boltzmann
        approximation, i.e. that the mass is large compared to the temperature. In this limit, even fewer terms of the
        expansion need to be kept. Doing so boosts performance.

        We work in the grand canonical ensemble, i.e. we consider P(V,T,mu_i/T). Hence in this class, derivatives w.r.t.
        one of those N_chemical_potentials + 2 variables assume all others are held fixed. """

    def __init__(self, Mass, g, w, B, S, Q, C=None, NMAX_light=21, NMAX_heavy=2):
        HRGbase.__init__(self, Mass, g, w, B, S, Q, C)
        self.NMAX_light = NMAX_light
        self.NMAX_heavy = NMAX_heavy

    def __repr__(self) -> str:
        return "HRG"

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
        """ (m/T)^2 g eta^(n+1) / 2pi^2 n^2 """
        return (self.Mass[k]/T)**2 * self.g[k] * self.w[k]**(n+1) / (2*np.pi**2*n**2)

    def P_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        P = 0.
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                P += self.factor(k, n, T) \
                        * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                        * kn(2, (n * self.Mass[k] / T))
        return P

    def E_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        eps = 0.
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                x = self.Mass[k]*n/T
                eps += self.factor(k, n, T) \
                        * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                        * (kn(2, x) * 3 + kn(1, x) * x)
        return eps

    def S_div_T3(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ s = e + p - mu_i n_i """
        NB = self.gen_chi(T, B_order=1, S_order=0, Q_order=0, C_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        NS = self.gen_chi(T, B_order=0, S_order=1, Q_order=0, C_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        NQ = self.gen_chi(T, B_order=0, S_order=0, Q_order=1, C_order=0, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        NC = self.gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=1, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        muxN_div_T = NB*muB_div_T + NS*muS_div_T + NQ*muQ_div_T + NC*muC_div_T
        return self.E_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T) \
                + self.P_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T,muQ_div_T=muQ_div_T, muC_div_T=muC_div_T) - muxN_div_T

    def ddT_E_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ d(E/T^4)/dT at fixed mu/T """
        eps = 0.
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                m = self.Mass[k]
                x = m*n/T
                eps += self.factor(k, n, T)*m*n \
                        * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                        * (kn(0, x)*n*m + kn(1, x)*T)/T**3
        return eps

    def ddT_P_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ d(P/T^4)/dT at fixed mu/T """
        P = 0.
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                m = self.Mass[k]
                x = m*n/T
                P += self.factor(k, n, T) \
                        * self.z(k, muB_div_T=muB_div_T,muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                        * m*n*kn(1, x) / T**2
        return P

    def CV_div_T3_mu0(self, T):
        """ C_V/T^3 at mu=0. """
        return 4*self.E_div_T4(T, 0, 0, 0, 0) + T*self.ddT_E_div_T4(T, 0, 0, 0, 0)

    def ddT_S_div_T3(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ d(s/T^3)/dT at fixed mu/T """
        ddTNB = self.ddT_gen_chi(T, B_order=1, S_order=0, Q_order=0, C_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        ddTNS = self.ddT_gen_chi(T, B_order=0, S_order=1, Q_order=0, C_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        ddTNQ = self.ddT_gen_chi(T, B_order=0, S_order=0, Q_order=1, C_order=0,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        ddTNC = self.ddT_gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=1,
                                 muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        muxddTN_div_T = ddTNB*muB_div_T + ddTNS*muS_div_T + ddTNQ*muQ_div_T + ddTNC*muC_div_T
        return self.ddT_E_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T) \
                + self.ddT_P_div_T4(T, muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T) - muxddTN_div_T

    def gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0.):
        """ chi_BQSC """
        chi = 0.0
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                zn_Kn = self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n * kn(2, (n*self.Mass[k]/T))
                chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order * (self.Q[k]*n)**Q_order * (self.C[k]*n)**C_order * self.factor(k, n, T) * zn_Kn
        return chi

    def ddT_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0.):
        """ d(chi_BQSC)/dT at fixed mu/T """
        chi = 0.0
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                m = self.Mass[k]
                x = m*n/T
                chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order * (self.Q[k]*n)**Q_order * (self.C[k]*n)**C_order \
                                              * self.factor(k, n, T) \
                                              * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                              * m*n*kn(1, x) / T**2
        return chi

    def d2dT2_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0.):
        """ d^2(chi_BQSC)/dT^2 at fixed mu/T """
        chi = 0.0
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                m = self.Mass[k]
                x = m*n/T
                chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order * (self.Q[k]*n)**Q_order * (self.C[k]*n)**C_order \
                                              * self.factor(k, n, T) * m*n \
                                              * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                              * (m*n*kn(0, x)/T**4 - 3*kn(1, x)/T**3)
        return chi

    def gen_ddmuh_E_div_T4(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ Arbitrary mu/T derivatives of E/T^4 at fixed T """
        eps = 0.
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                x = self.Mass[k]*n/T
                eps += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order * (self.Q[k]*n)**Q_order * (self.C[k]*n)**C_order \
                                              * self.factor(k, n, T) \
                                              * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                              * (kn(2, x)*3 + kn(1, x)*x)
        return eps

    def gen_ddmuh_P_div_T4(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ Arbitrary mu/T derivatives of P/T^4 at fixed T """
        P = 0.
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                P += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order * (self.Q[k]*n)**Q_order * (self.C[k]*n)**C_order \
                                            * self.factor(k, n, T) \
                                            * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                            * kn(2, (n*self.Mass[k]/T))
        return P

    def gen_ddmuh_S_div_T3(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ s = e + p - mu_i n_i """
        NB = self.gen_chi(T, B_order=1, S_order=0, Q_order=0, C_order=0, muB_div_T=muB_div_T,
                          muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        NS = self.gen_chi(T, B_order=0, S_order=1, Q_order=0, C_order=0, muB_div_T=muB_div_T,
                          muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        NQ = self.gen_chi(T, B_order=0, S_order=0, Q_order=1, C_order=0, muB_div_T=muB_div_T,
                          muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        NC = self.gen_chi(T, B_order=0, S_order=0, Q_order=0, C_order=1, muB_div_T=muB_div_T,
                          muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        dmuhNB = self.gen_chi(T, B_order=1+B_order, S_order=S_order, Q_order=Q_order, C_order=C_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        dmuhNS = self.gen_chi(T, B_order=B_order, S_order=1+S_order, Q_order=Q_order, C_order=C_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        dmuhNQ = self.gen_chi(T, B_order=B_order, S_order=S_order, Q_order=1+Q_order, C_order=C_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        dmuhNC = self.gen_chi(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=1+C_order,
                              muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)
        return self.gen_ddmuh_E_div_T4(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=C_order, 
                                       muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T) \
                + self.gen_ddmuh_P_div_T4(T, B_order=B_order, S_order=S_order, Q_order=Q_order, C_order=C_order, 
                                          muB_div_T=muB_div_T, muS_div_T=muS_div_T, muQ_div_T=muQ_div_T, muC_div_T=muC_div_T) \
                - muB_div_T*dmuhNB - muQ_div_T*dmuhNQ - muS_div_T*dmuhNS - muC_div_T*dmuhNC \
                - NB*dmuh(B_order, muB_div_T) - NQ*dmuh(Q_order, muQ_div_T) - NS*dmuh(S_order, muS_div_T) - NC*dmuh(C_order, muC_div_T)

    def gen_chi_RMS(self, T, Nt, B_order=0, S_order=0, Q_order=0, C_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0.):
        # rms_mass[0] is for pions and rms_mass[1] is for kaons
        rms_mass = RMS_mass(Nt, T)
        chi = 0.0
        for k in range(len(self.Mass)):
            if 140 >= self.Mass[k] >= 130:
                for n in range(1, 20):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order \
                                                  * (self.Q[k]*n)**Q_order \
                                                  * (self.C[k]*n)**C_order \
                                                  * self.w[k]**(n+1) * self.g[k] * (rms_mass[0]/T)**2 \
                                                  * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                                  * kn(2, (n*rms_mass[0]/T)) / (np.pi*n)**2 / 2
            elif 500 >= self.Mass[k] >= 490:
                for n in range(1, 10):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order \
                                                  * (self.Q[k]*n)**Q_order \
                                                  * (self.C[k]*n)**C_order \
                                                  * self.w[k]**(n+1) * self.g[k] * (rms_mass[1]/T)**2 \
                                                  * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                                  * kn(2, (n*rms_mass[1]/T)) / (np.pi*n)**2 / 2
            else:
                for n in range(1, 2):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order \
                                                  * (self.Q[k]*n)**Q_order \
                                                  * (self.C[k]*n)**C_order \
                                                  * self.factor(k, n, T) \
                                                  * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                                                  * kn(2, (n*self.Mass[k]/T))
        return chi

    def genChiFlavor(self, T, u_order=0, d_order=0, s_order=0, c_order=0, muB_div_T=0., muQ_div_T=0., muS_div_T=0., muC_div_T=0.):
        """ chi_udsc """
        logger.warn('HRGFlavor basis is not yet tested for finite mu...')
        chi = 0.0
        for k in range(len(self.Mass)):
            for n in range(1, self.Nmax(k)):
                chi += ((self.B[k]+self.Q[k]-self.C[k])*n)**u_order \
                    * ((2.0*self.B[k]-self.Q[k]+self.S[k])*n)**d_order \
                    * (-self.S[k]*n)**s_order \
                    * (self.C[k]*n)**c_order \
                    * self.factor(k, n, T) \
                    * self.z(k, muB_div_T=muB_div_T, muQ_div_T=muQ_div_T, muS_div_T=muS_div_T, muC_div_T=muC_div_T)**n \
                    * kn(2, (n*self.Mass[k]/T))
        return chi
# TODO: gen_chi can be implemented using multi-dimensional Leibniz rule


class HRGexact(HRGbase):
# TODO: maybe just delete this eventually, it seems to be crap

    """ HRG implemented through numerical integration. """

    def __init__(self, Mass, g, w, B, S, Q, C=None):
        logger.warn('ExactHRG numerical integration is not yet reliable...')
        HRGbase.__init__(self, Mass, g, w, B, S, Q, C)

    def __repr__(self) -> str:
        return "HRGexact"

    def P_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T = envector(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)

        def int_wrapper(Tvec, muBvec, muSvec, muQvec, muCvec):
            P = 0.
            for k in range(len(self.Mass)):
                wz = self.w[k] * self.z(k, muB_div_T=muBvec,muQ_div_T=muQvec, muS_div_T=muSvec, muC_div_T=muCvec)
                def integrand(E):
                    exp_E_div_T = np.exp(-E/Tvec)
                    return -(wz/(3*Tvec)) * (E**2-self.Mass[k]**2)**(3/2) * exp_E_div_T / (1-wz*exp_E_div_T)
                P -= self.w[k] * self.g[k] * integrateFunction(integrand, self.Mass[k], np.inf, method='quad') / (2*np.pi**2 * Tvec**3)
            return P
        int_vec = np.vectorize(int_wrapper)
        return unvector(np.asarray(int_vec(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)))

    def E_div_T4(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T = envector(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)

        def int_wrapper(Tvec, muBvec, muSvec, muQvec, muCvec):
            eps = 0.
            for k in range(len(self.Mass)):
                wz = self.w[k] * self.z(k, muB_div_T=muBvec,
                                        muQ_div_T=muQvec, muS_div_T=muSvec, muC_div_T=muCvec)

                def integrand(E):
                    exp_E_div_T = np.exp(-E/Tvec)
                    return wz * E**2 * (E**2-self.Mass[k]**2)**(1/2) * exp_E_div_T / (1-wz*exp_E_div_T)
                eps += self.w[k] * self.g[k] * integrateFunction(integrand, self.Mass[k], np.inf, method='quad')
            eps /= (2*np.pi**2*Tvec**4)
            return eps
        int_vec = np.vectorize(int_wrapper)
        return unvector(np.asarray(int_vec(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)))

    def number_density(self, T, charge='B', muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T = envector(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)
        if charge == 'B':
            X = self.B
        elif charge == 'Q':
            X = self.Q
        elif charge == 'S':
            X = self.S
        elif charge == 'C':
            X = self.C
        else:
            logger.TBError('Unrecognized charge', charge)

        def int_wrapper(Tvec, muBvec, muSvec, muQvec, muCvec):
            NX = 0.
            for k in range(len(self.Mass)):
                wz = self.w[k] * self.z(k, muB_div_T=muBvec,muQ_div_T=muQvec, muS_div_T=muSvec, muC_div_T=muCvec)
                def integrand(E):
                    exp_E_div_T = np.exp(-E/Tvec)
                    return -wz * X[k] * E * (E**2-self.Mass[k]**2)**(1/2) * exp_E_div_T / (1 - wz*exp_E_div_T)
                NX -= self.w[k] * self.g[k] * integrateFunction(integrand, self.Mass[k], np.inf, method='quad')
            NX /= (2*np.pi**2*Tvec**3)
            return NX
        int_vec = np.vectorize(int_wrapper)
        return unvector(np.asarray(int_vec(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)))

    def S_div_T3(self, T, muB_div_T=0., muS_div_T=0., muQ_div_T=0., muC_div_T=0.):
        """ s = e + p - mu_i n_i """
        T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T = envector(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)

        def int_wrapper(Tvec, muBvec, muSvec, muQvec, muCvec):
            S = 0.
            for k in range(len(self.Mass)):
                wz = self.w[k] * self.z(k, muB_div_T=muBvec, muQ_div_T=muQvec, muS_div_T=muSvec, muC_div_T=muCvec)

                def integrand(E):
                    exp_E_div_T = np.exp(-E/Tvec)
                    P = (E**2-self.Mass[k]**2)/3
                    eps = E**2
                    return ( wz * (E**2-self.Mass[k]**2)**(1/2) \
                                * (P + eps - Tvec*E*(self.B[k]*muBvec + self.Q[k]*muQvec + self.S[k]*muSvec + self.C[k]*muCvec)) \
                                * exp_E_div_T / (1-wz*exp_E_div_T) )

                S += self.w[k] * self.g[k] * integrateFunction(integrand, self.Mass[k], np.inf, method='quad')
            S /= (2*np.pi**2*Tvec**4)
            return S
        int_vec = np.vectorize(int_wrapper)
        return unvector(np.asarray(int_vec(T, muB_div_T, muS_div_T, muQ_div_T, muC_div_T)))


class EVHRG(HRGbase):

    """ Excluded volume hadron resonance gas. b is excluded volume parameter. """

    def __init__(self, Mass, g, w, B, S, Q, C=None):
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
