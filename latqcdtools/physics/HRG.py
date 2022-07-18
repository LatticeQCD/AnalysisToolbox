#
# HRG.py
#
# J. Goswami, D. Clarke
#
# Collection of methods pertaining to hadron resonance gas calculations.
#

import numpy as np
import sympy as sy
from scipy.special import kn, lambertw
from sympy import Sum, symbols, Indexed, lambdify, LambertW, exp


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


class HRG:

    """ Hadron resonance gas. Mass=mass of the Hadron , g=spin degenerecy , w= fermi(-1)/bose(1) statistics.
        B, Q, S, and C are respectively the baryon number, electric charge, strangeness, and charm of each state. For
        more information please see, e.g. Physics Letters B 695 (2011) 136–142 or especially arXiv:2011.02812.

        Our pressure is given in terms of a Taylor series involving modified Bessel functions of the second kind, which
        needs to be truncated at some order. These functions get strongly suppressed when their argument is large.
        In our case, this argument is proportional to the mass. Hence we will sometimes take the Boltzmann
        approximation, i.e. that the mass is large compared to the temperature. In this limit, even fewer terms of the
        expansion need to be kept. Doing so boosts performance. """


    # For now keep parallelize false. This implementation does not speed anything up for some reason.
    def __init__(self, Mass, g, w, B, S, Q, C = None):
        # M = Mass of the hadron
        # Q = charge of the hadron
        # B = Baryon number of the hadron [B=1,-1,0,0 for baryon,anti-baryon,meson,anti-mesons]
        # S = Strangeness number of the hadron
        # C = Charm number of the hadron
        # g = degenracy of the hadron state
        # w = spin statistics of hadron (eta)
        self.Mass = Mass
        self.g = g
        self.w = w
        self.B = B
        self.Q = Q
        self.S = S
        # If we don't get a charm array, initialize to array of zeroes.
        if C is None:
            self.C = np.zeros(len(Mass))
        else:
            self.C = C


    def __repr__(self):
        return "hrg"


    # n represents the nth order of a Taylor expansion of a logarithm.
    # k represents the kth state that appears in the table of resonances.
    def factor(self, k, n, T):
        # m^2 g eta^(n+1) T^2 / 2pi^2 n^2
        return self.Mass[k]**2 * self.g[k] * T**2 * self.w[k]**(n+1) / (2*np.pi**2*n**2)


    def muN(self, k, mu_B, mu_Q, mu_S, mu_C):
        # mu_i * N_i
        return self.B[k]*mu_B + self.Q[k]*mu_Q + self.S[k]*mu_S + self.C[k]*mu_C


    def z(self, T, k, mu_B, mu_Q, mu_S, mu_C):
        # e^(mu_i*N_i/T)
        return np.exp( self.muN(k, mu_B, mu_Q, mu_S, mu_C)/T )


    def pressure(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        P = 0.0
        for k in range(len(self.Mass)):
            for n in range(1, 20): # Keep only first 20 terms of the series.
                P += self.factor(k, n, T) * self.z(T, k, mu_B, mu_Q, mu_S, mu_C)**n * kn(2,(n*self.Mass[k]/T))
        return P


    def P_div_T4(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        return self.pressure(T, mu_B, mu_S, mu_Q, mu_C)/T**4


    def energy_density(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        eps = 0.
        for k in range(len(self.Mass)):
            for n in range(1,20):
                x = self.Mass[k]*n/T
                eps += self.factor(k,n,T) * self.z(T,k,mu_B,mu_Q,mu_S,mu_C)**n \
                       * ( kn(2,x) * 3 + kn(1,x)*x )
        return eps


    def E_div_T4(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        return self.energy_density(T, mu_B, mu_S, mu_Q, mu_C)/T**4


    # all unitless: s = e + p - mu_i n_i
    def entropy(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        muxN=0.
        for k in range(len(self.Mass)):
            muxN += self.muN(k, mu_B, mu_Q, mu_S, mu_C)
        return ( self.energy_density(T,mu_B,mu_S,mu_Q,mu_C) + self.pressure(T,mu_B,mu_S,mu_Q,mu_C) - muxN )/T


    def S_div_T3(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        return self.entropy(T,mu_B,mu_S,mu_Q,mu_C)/T**3


    def ddT_E_div_T4(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        eps = 0.
        for k in range(len(self.Mass)):
            for n in range(1,20):
                m    = self.Mass[k]
                x    = m*n/T
                eps += self.factor(k,n,T)*m*n * self.z(T,k,mu_B,mu_Q,mu_S,mu_C)**n * ( kn(0,x)*n*m + kn(1,x)*T )/T**7
        return eps


    def ddT_P_div_T4(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        P = 0.
        for k in range(len(self.Mass)):
            muxN = self.muN(k, mu_B, mu_Q, mu_S, mu_C)
            for n in range(1,20):
                m    = self.Mass[k]
                x    = m*n/T
                P += self.factor(k, n, T) * n * self.z(T, k, mu_B, mu_Q, mu_S, mu_C)**n * ( kn(1,x) * m - kn(2,x) * muxN )/T**2
        return P


    def CV(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        return ( 4*self.E_div_T4(T, mu_B, mu_S, mu_Q, mu_C) + T*self.ddT_E_div_T4(T, mu_B, mu_S, mu_Q, mu_C) ) * T**3


    def CV_div_T3(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        return self.CV(T, mu_B, mu_S, mu_Q, mu_C)/T**3


    def gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, mu_B=0., mu_Q=0., mu_S=0., mu_C=0.):
        chi = 0.0
        for k in range(len(self.Mass)):
            Nterms = 20
            if self.B[k] != 0:
                Nterms = 2
            for n in range(1, Nterms):
                chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order \
                                              * (self.Q[k]*n)**Q_order \
                                              * (self.C[k]*n)**C_order \
                                              * self.factor(k, n, T) * self.z(T, k, mu_B, mu_Q, mu_S, mu_C)**n \
                                              * kn(2,(n*self.Mass[k]/T))
        return chi/T**4


    def gen_chi_RMS(self, T, Nt, B_order=0, S_order=0, Q_order=0, C_order=0, mu_B=0., mu_Q=0., mu_S=0., mu_C=0.):
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
                                                  * self.z(T, k, mu_B, mu_Q, mu_S, mu_C)**n \
                                                  * kn(2, (n*rms_mass[0]/T)) / (np.pi*n)**2 / 2
            elif 500 >= self.Mass[k] >= 490:
                for n in range(1, 10):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order \
                                                  * (self.Q[k]*n)**Q_order \
                                                  * (self.C[k]*n)**C_order \
                                                  * self.w[k]**(n+1) * self.g[k] * (rms_mass[1]/T)**2 \
                                                  * self.z(T, k, mu_B, mu_Q, mu_S, mu_C)**n \
                                                  * kn(2, (n*rms_mass[1]/T)) / (np.pi*n)**2 / 2
            else:
                for n in range(1, 2):
                    chi += (self.B[k]*n)**B_order * (self.S[k]*n)**S_order \
                                                  * (self.Q[k]*n)**Q_order \
                                                  * (self.C[k]*n)**C_order \
                                                  * self.factor(k, n, T) * self.z(T, k, mu_B, mu_Q, mu_S, mu_C)**n \
                                                  * kn(2,(n*self.Mass[k]/T))/T**4
        return chi


# TODO: the __init__ can be inherited from the HRG class
class EV_HRG:

    """ Excluded volume hadron resonance gas. Mass=mass of the Hadron , g=spin degenerecy , w= fermi(-1)/bose(1) statistics. """

    def __init__(self, Mass, g, w, B, S, Q):
        self.Mass = Mass
        self.g = g
        self.w = w
        self.B = B
        self.S = Q
        self.Q = S

    def __repr__(self):
        return "evhrg"

    def Pid(self, m, g, T):
        return g*(m/T)**2 * kn(2, (m/T)) / np.pi**2 / 2

    def baryon_pressure(self, T, b, Bi, mu_B=0.0, mu_Q=0.0, mu_S=0.0):

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

        mB, mQ, mS, i, ci, bi, qi, si, temp = symbols('mB, mQ , mS, i, ci, bi, qi, si, temp')

        F_pressure = (1 / (b*(T/197.3)**3)) * LambertW(
            (b*(T/197.3)**3) * Sum(Indexed('ci', i) * exp(mB * Indexed('bi', i) + mQ * Indexed('qi', i) + mS * Indexed('si', i)),
                                   (i, 0, len(P) - 1)))
        f = lambdify((mB, mQ, mS, ci, bi, qi, si, temp), F_pressure,
                     modules=['scipy', {'LambertW': lambertw}])

        pressure = f(mu_B/T, mu_Q/T, mu_S/T, P, X_baryon, X_charge, X_strange, T).real

        return pressure

    def gen_chi(self, T, b, Bi, B_order=0, Q_order=0, S_order=0, mu_B = 0.0, mu_Q = 0.0, mu_S = 0.0):

        baryon_mass = self.Mass[np.where(self.B == Bi)]
        g_baryon    = self.g[np.where(self.B    == Bi)]
        X_baryon    = self.B[np.where(self.B    == Bi)]
        X_charge    = self.Q[np.where(self.B    == Bi)]
        X_strange   = self.S[np.where(self.B    == Bi)]

        # Bi=1 Baryon pressure , Bi=-1 for anti baryon pressure
        P = []
        for k in range(len(baryon_mass)):
            P.append(self.Pid(baryon_mass[k], g_baryon[k], T))

        P = np.array(P)

        # ci = ideal gas pressure of individual particles
        mB, mQ, mS, i, ci, bi, qi, si, temp, be = symbols('mB, mQ , mS, i, ci, bi, qi, si, temp, be')

        F_pressure = (1 / be) * LambertW(
            be * Sum(
                Indexed('ci', i) * exp(mB * Indexed('bi', i) + mQ *
                                       Indexed('qi', i) + mS * Indexed('si', i)),
                (i, 0, len(P) - 1)))
        expr_chi = sy.diff(F_pressure, mB, B_order, mQ, Q_order, mS, S_order)

        f = lambdify((mB, mQ, mS, ci, bi, qi, si, temp, be),
                     expr_chi, modules=['scipy', {'LambertW': lambertw}])

        chi_num = f(mu_B/T, mu_Q/T, mu_S/T, P, X_baryon,
                    X_charge, X_strange, T, b*(T/197.3)**3).real

        return chi_num
