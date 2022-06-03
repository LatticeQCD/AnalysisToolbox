# 
# HRG.py                                                               
# 
# J. Goswami 
#
# Collection of methods pertaining to hadron resonance gas calculations.
#


import numpy as np
import sympy as sy
from scipy.special import kn, lambertw
from sympy import Sum, symbols, Indexed, lambdify, LambertW, exp


def RMS_mass(Nt,T):

    a  = 6924.46
    b  = 31881.4
    c  = -1.02357e+06
    mkaon=493.68
    mpion=140.0

    x=1.0/Nt/T*197.3
    pion_rms_mass = mpion + a*x**2+b*x**4+c*x**6

    del2= pion_rms_mass**2 - mpion**2
    kaon_rms_mass= np.sqrt(mkaon**2+del2)

    return pion_rms_mass,kaon_rms_mass


class HRG:

    """ Hadron resonance gas. Mass=mass of the Hadron , g=spin degenerecy , w= fermi(-1)/bose(1) statistics.
        B = baryon number of HRG state, Q = charge HRG state, S = strangeness HRG state. """

    def __init__(self, Mass, g, w, B=0.0, S=0.0, Q=0.0):
        self.Mass = Mass
        self.g = g
        self.w = w
        self.B = B
        self.S = S
        self.Q = Q

    def __repr__(self):
        return "hrg"

    def ln_Z(self, k, N, T):
        return self.w[k]**(N + 1)*self.g[k]*(self.Mass[k]/T)**2*kn(2, (N * self.Mass[k]/T))/(np.pi*N)**2/2

    def exp(self,N,T,k,mu_B,mu_Q,mu_S):
         return np.exp(N * (self.B[k] * mu_B + self.Q[k] * mu_Q + self.S[k] * mu_S) / T)

    def pressure(self, T, mu_B=0., mu_S=0., mu_Q=0.):
        P = 0.0
        for k in range(len(self.Mass)):
            for N in range(1, 20):
                if N * self.Mass[k] > 2500:
                    y = 0.0
                else:
                    y = self.ln_Z(k, N, T)*np.exp( N* (self.B[k]*mu_B + self.Q[k]*mu_Q + self.S[k]*mu_S)/ T )
                P += y
        return P

    def gen_chi(self, T, B_order = 0.0, S_order = 0.0, Q_order = 0.0 ,mu_B = 0.0, mu_Q = 0.0, mu_S = 0.0):
        chi = 0.0
        for k in range(len(self.Mass)):
            if self.B[k]==0:
                for N in range(1,20):
                    chi += (self.B[k] * N) ** B_order * (self.S[k] * N) ** S_order * (self.Q[k] * N) ** Q_order * self.ln_Z(k, N, T) * np.exp(N * (self.B[k] * mu_B + self.Q[k] * mu_Q + self.S[k] * mu_S) / T)
            else:
                for N in range(1,2):
                    chi += (self.B[k] * N) ** B_order * (self.S[k] * N) ** S_order * (self.Q[k] * N) ** Q_order * self.ln_Z(k, N, T) * np.exp(N * (self.B[k] * mu_B + self.Q[k] * mu_Q + self.S[k] * mu_S) / T)
        return chi

    def gen_chi_RMS(self, T, Nt,B_order=0.0, S_order=0.0, Q_order=0.0 , mu_B = 0.0, mu_Q = 0.0, mu_S = 0.0 ):
        rms_mass=RMS_mass(Nt,T) #rms_mass[0] is for pions and rms_mass[1] is for kaons
        chi = 0.0
        for k in range(len(self.Mass)):
            if 140 >= self.Mass[k] >= 130:
                for N in range(1,20):
                    chi += (self.B[k] * N) ** B_order * (self.S[k] * N) ** S_order * (self.Q[k] * N) ** Q_order * self.w[k] ** (N + 1) * self.g[k] * (rms_mass[0] / T) ** 2 * self.exp(N , T , k , mu_B , mu_Q , mu_S) * kn(2, (
                N * rms_mass[0] / T)) / (np.pi * N) ** 2 / 2
            elif 500 >= self.Mass[k] >= 490:
                for N in range(1,10):
                    chi += (self.B[k] * N) ** B_order * (self.S[k] * N) ** S_order * (self.Q[k] * N) ** Q_order * self.w[k] ** (N + 1) * self.g[k] * (rms_mass[1] / T) ** 2 * self.exp(N , T , k , mu_B , mu_Q , mu_S) * kn(2, (
                N * rms_mass[1] / T)) / (np.pi * N) ** 2 / 2
            else:
                for N in range(1,2):
                    chi += (self.B[k] * N) ** B_order * (self.S[k] * N) ** S_order * (self.Q[k] * N) ** Q_order * self.ln_Z(k, N, T) * self.exp(N , T , k , mu_B , mu_Q , mu_S)
        return chi


class EV_HRG:

    """ Excluded volume hadron resonance gas. Mass=mass of the Hadron , g=spin degenerecy , w= fermi(-1)/bose(1) statistics. """

    def __init__(self, Mass, g, w, B=0.0, S=0.0, Q=0.0):
        self.Mass = Mass
        self.g = g
        self.w = w
        self.B = B
        self.S = S
        self.Q = Q

    def __repr__(self):
        return "evhrg"

    def Pid(self, m, g, T):
        return  g*(m/T)**2  * kn(2, (m/T)) / np.pi**2 / 2

    def exp(self,N,T,k,mu_B,mu_Q,mu_S):
         return np.exp(N * (self.B[k] * mu_B + self.Q[k] * mu_Q + self.S[k] * mu_S) / T)

    def baryon_pressure(self, T, b, Bi, mu_B=0.0, mu_Q=0.0, mu_S=0.0):

        baryon_mass = self.Mass[np.where(self.B == Bi)]
        g_baryon    = self.g[np.where(self.B == Bi)]
        X_baryon    = self.B[np.where(self.B == Bi)]
        X_charge    = self.Q[np.where(self.B == Bi)]
        X_strange   = self.S[np.where(self.B == Bi)]

        # Bi=1 Baryon pressure , Bi=-1 for anti baryon pressure
        P=[]
        for k in range(len(baryon_mass)):
            P.append(self.Pid(baryon_mass[k],g_baryon[k], T))

        P= np.array(P)

        mB, mQ, mS, i, ci, bi, qi, si, temp = symbols('mB, mQ , mS, i, ci, bi, qi, si,temp')

        F_pressure = (1 / (b*(T/197.3)**3)) * LambertW(
        (b*(T/197.3)**3) * Sum(Indexed('ci', i) * exp(mB * Indexed('bi', i) + mQ * Indexed('qi', i) + mS * Indexed('si', i)),
                         (i, 0, len(P) - 1)))
        f = lambdify((mB, mQ, mS, ci, bi, qi, si, temp), F_pressure, modules=['scipy', {'LambertW': lambertw}])

        pressure = f(mu_B/T, mu_Q/T, mu_S/T, P, X_baryon, X_charge, X_strange, T).real

        return pressure

    def gen_chi(self, T, b, Bi, B_order=0, Q_order=0, S_order=0,mu_B=0.0, mu_Q=0.0, mu_S=0.0):

        baryon_mass = self.Mass[np.where(self.B == Bi)]
        g_baryon    = self.g[np.where(self.B == Bi)]
        X_baryon    = self.B[np.where(self.B == Bi)]
        X_charge    = self.Q[np.where(self.B == Bi)]
        X_strange   = self.S[np.where(self.B == Bi)]

        # Bi=1 Baryon pressure , Bi=-1 for anti baryon pressure
        P = []
        for k in range(len(baryon_mass)):
            P.append(self.Pid(baryon_mass[k], g_baryon[k], T))

        P = np.array(P)

        # ci = ideal gas pressure of individual particles
        mB, mQ, mS, i, ci, bi, qi, si, temp,be = symbols('mB, mQ , mS, i, ci, bi, qi, si, temp, be')

        F_pressure = (1 / be) * LambertW(
             be * Sum(
                Indexed('ci', i) * exp(mB * Indexed('bi', i) + mQ * Indexed('qi', i) + mS * Indexed('si', i)),
                (i, 0, len(P) - 1)))
        expr_chi = sy.diff(F_pressure, mB, B_order, mQ, Q_order, mS, S_order)

        f = lambdify((mB, mQ, mS, ci, bi, qi, si, temp, be), expr_chi, modules=['scipy', {'LambertW': lambertw}])

        chi_num = f(mu_B/T, mu_Q/T, mu_S/T, P, X_baryon, X_charge, X_strange, T, b*(T/197.3)**3).real

        return chi_num


