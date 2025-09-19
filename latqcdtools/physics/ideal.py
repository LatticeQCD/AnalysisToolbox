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
import latqcdtools.base.logger as logger

class idealGas:
    """ 
    Class for ideal, massless gas of quarks, zeroth order in the coupling. Fermion chemical
    potentials are added in the order u-->d-->s-->c.

    Args:
        Nc (int): Number of colors. 
        Nf (int): Number of fermion flavors. 
    """

    def __init__(self,Nc,Nf): 
        checkType('int',Nc=Nc)
        checkType('int',Nf=Nf)
        if Nf >= 4:
            logger.TBRaise('Only B, Q, S, C explicitly coded in.')
        if Nc < 1:
            logger.TBRaise('Must have positive number of colors.')
        self.T, self.muB, self.muS, self.muQ, self.muC = sympy.symbols('T muB muS muQ muC') 
        self.Nf   = Nf
        self.Nc   = Nc
        self.muu  = (self.muB + 2*self.muQ)/3
        self.mud  = (self.muB -   self.muQ)/3
        self.mus  = (self.muB -   self.muQ)/3 - self.muS
        self.muc  = (self.muB + 2*self.muQ)/3 + self.muC
        # eq (8.42) and (8.43) in Kapusta and Gale's "Finite temperature field theory: 
        # Principles and Applications" (2006)
        self.Psym = (np.pi**2/45)*(self.Nc**2-1)*self.T**4 
        if Nf > 0:  
            self.Psym += self.Nc*( 7*np.pi**2*self.T**4/180 + self.muu**2*self.T**2/6 + self.muu**4/(12*np.pi**2) )
        if Nf > 1:  
            self.Psym += self.Nc*( 7*np.pi**2*self.T**4/180 + self.mud**2*self.T**2/6 + self.mud**4/(12*np.pi**2) )
        if Nf > 2:  
            self.Psym += self.Nc*( 7*np.pi**2*self.T**4/180 + self.mus**2*self.T**2/6 + self.mus**4/(12*np.pi**2) )
        if Nf > 3:  
            self.Psym += self.Nc*( 7*np.pi**2*self.T**4/180 + self.muc**2*self.T**2/6 + self.muc**4/(12*np.pi**2) )

    def __repr__(self) -> str:
        return "idealGas"

    def P(self, T, muB=0., muS=0., muQ=0., muC=0.):
        """ 
        Unitful pressure. 
        """
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS, self.muC: muC} 
        return np.float128(self.Psym.subs(values).evalf())

    def gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB=0., muQ=0., muS=0., muC=0.):
        """
        Unitful cumulant.
        """
        chi = sympy.diff(self.Psym, self.muB, B_order)
        chi = sympy.diff(chi      , self.muQ, Q_order)
        chi = sympy.diff(chi      , self.muS, S_order)
        chi = sympy.diff(chi      , self.muC, C_order)
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS, self.muC: muC} 
        return np.float128(chi.subs(values).evalf())

    def S(self, T, muB=0., muS=0., muQ=0., muC=0.):
        """ 
        Unitful entropy. 
        """
        entropy = sympy.diff(self.Psym, self.T, 1)
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS, self.muC: muC} 
        return np.float128(entropy.subs(values).evalf())

#    def E(self, T, muB=0., muS=0., muQ=0., muC=0.):
#        """
#        Unitful energy density.
#        """
#        epsilon = self.T**2*sympy.diff(self.Psym/self.T, self.T, 1)
#        values  = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS, self.muC: muC} 
#        return float(epsilon.subs(values).evalf())

    def ddT_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB=0., muQ=0., muS=0., muC=0.):
        chi     = sympy.diff(self.Psym, self.muB, B_order)
        chi     = sympy.diff(chi      , self.muQ, Q_order)
        chi     = sympy.diff(chi      , self.muS, S_order)
        chi     = sympy.diff(chi      , self.muC, C_order)
        ddT_chi = sympy.diff(chi      , self.T  , 1      )
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS, self.muC: muC} 
        return np.float128(ddT_chi.subs(values).evalf())

    def d2dT2_gen_chi(self, T, B_order=0, S_order=0, Q_order=0, C_order=0, muB=0., muQ=0., muS=0., muC=0.):
        chi       = sympy.diff(self.Psym, self.muB, B_order)
        chi       = sympy.diff(chi      , self.muQ, Q_order)
        chi       = sympy.diff(chi      , self.muS, S_order)
        chi       = sympy.diff(chi      , self.muC, C_order)
        d2dT2_chi = sympy.diff(chi      , self.T  , 2      )
        values = {self.T: T, self.muB: muB, self.muQ: muQ, self.muS: muS, self.muC: muC} 
        return np.float128(d2dT2_chi.subs(values).evalf())

