#
# ideal.py 
#
# D. Clarke 
#
# Methods for ideal gas calculations. 
#


import numpy as np
import sympy
from latqcdtools.base.check import checkType


class idealGas:

    def __init__(self,nf): 
        """ The idealGas class. Analytic derivatives are implemented with sympy.

        Args:
            nf (int): Number of fermion flavors. 
        """
        checkType(nf,int)
        self.T, self.muB, self.muS, self.muQ = sympy.symbols('T muB muS muQ') 
        self.nf   = nf
        self.c0   = 8*np.pi**2/45*(1 + 21*self.nf/32)
        self.muu  = (self.muB + 2*self.muQ)/3
        self.mud  = (self.muB -   self.muQ)/3
        self.mus  = (self.muB -   self.muQ)/3 - self.muS
        self.Psym = self.T**4*self.c0 + self.T**2*(self.muu**2 + self.mud**2 + self.mus**2)/2 \
                    + (self.muu**4 + self.mud**4 + self.mus**4)/(4*np.pi**2)

    def __repr__(self) -> str:
        return "idealGas"

    def P(self, T, muB=0., muS=0., muQ=0., muC=0.):
        """ Unitful pressure. """
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS} 
        return float(self.Psym.subs(values).evalf())

    def gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB=0., muQ=0., muS=0., muC=0.):
        chi = sympy.diff(self.Psym, self.muB, B_order)
        chi = sympy.diff(chi      , self.muQ, Q_order)
        chi = sympy.diff(chi      , self.muS, S_order)
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS} 
        return float(chi.subs(values).evalf())

    def S(self, T, muB=0., muS=0., muQ=0., muC=0.):
        """ Unitful entropy. """
        entropy = sympy.diff(self.Psym, self.T, 1)
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS} 
        return float(entropy.subs(values).evalf())

    def ddT_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB=0., muQ=0., muS=0., muC=0.):
        chi     = sympy.diff(self.Psym, self.muB, B_order)
        chi     = sympy.diff(chi      , self.muQ, Q_order)
        chi     = sympy.diff(chi      , self.muS, S_order)
        ddT_chi = sympy.diff(chi      , self.T  , 1      )
        values  = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS} 
        return float(ddT_chi.subs(values).evalf())

    def d2dT2_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB=0., muQ=0., muS=0., muC=0.):
        chi       = sympy.diff(self.Psym, self.muB, B_order)
        chi       = sympy.diff(chi      , self.muQ, Q_order)
        chi       = sympy.diff(chi      , self.muS, S_order)
        d2dT2_chi = sympy.diff(chi      , self.T  , 2      )
        values    = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS} 
        return float(d2dT2_chi.subs(values).evalf())

