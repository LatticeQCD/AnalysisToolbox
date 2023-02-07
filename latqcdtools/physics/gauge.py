# 
# gauge.py                                                               
# 
# D. Clarke
# 
# A gauge field class.
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.math.SU3 import SU3

class gaugeField:

    """ A class for SU3 gauge fields. """

    def __init__(self, Ns=None, Nt=None):
        """ Initialize field to unit matrix at every site. """
        self.Ns = Ns       # Spatial extension
        self.Nt = Nt       # Euclidean time extension
        self.field = np.zeros(shape=(self.Nt,self.Ns,self.Ns,self.Ns,4,3,3),dtype=complex)
        if (Ns is None) or (Nt is None):
           logger.TBError("gaugeField objects must be initialized with a size!")
        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(1,4):
                            U = SU3()
                            self.field[t,z,y,x,mu] = U


    def respectBCs(self,x,y,z,t):
        """ Assume periodic BCs. This also seems to work for negative numbers, e.g. -1 % 8 = 7. """
        return x%self.Ns, y%self.Ns, z%self.Ns, t%self.Nt


    def getLink(self,x,y,z,t,mu):
        """ Link accessor. """
        return self.field[t][z][y][x][mu]


    def setLink(self, other, x, y, z, t, mu):
        """ Link setter. """
        self.field[t][z][y][x][mu] = other


    def stepUp(self,x,y,z,t,mu):
        X, Y, Z, T = x, y, z, t
        if mu==0:
            X += 1
        elif mu==1:
            Y += 1
        elif mu==2:
            Z += 1
        elif mu==3:
            T += 1
        else:
            logger.TBError('mu =',mu,'violates 4d lattice.')
        return self.respectBCs(X, Y, Z, T)


    def stepDown(self,x,y,z,t,mu):
        X, Y, Z, T = x, y, z, t
        if mu==0:
            X -= 1
        elif mu==1:
            Y -= 1
        elif mu==2:
            Z -= 1
        elif mu==3:
            T -= 1
        else:
            logger.TBError('mu =',mu,'violates 4d lattice.')
        return self.respectBCs(X, Y, Z, T)


    def getLocalPlaquette(self,x,y,z,t,mu,nu):
        """ Get Re tr U_{mu,nu}^[](x,y,z,t). """
        U_mu_x     = SU3()
        U_nu_xpamu = SU3()
        U_mu_xpanu = SU3()
        U_nu_x     = SU3()
        U_mu_x.setToMatrix( self.getLink(x,y,z,t,mu) )
        X, Y, Z, T = self.stepUp(x,y,z,t,mu)
        U_nu_xpamu.setToMatrix( self.getLink(X,Y,Z,T,nu) )
        X, Y, Z, T = self.stepUp(x,y,z,t,nu)
        U_mu_xpanu.setToMatrix( self.getLink(X,Y,Z,T,mu) )
        U_nu_x.setToMatrix( self.getLink(x,y,z,t,nu) )
        plaq = U_mu_x * U_nu_xpamu * U_mu_xpanu.dagger() * U_nu_x.dagger()
        return plaq.trace().real


    def getPlaquette(self):
        """ Calculate <Re tr U_{mu,nu}^[]>. """
        plaq = 0.
        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for nu in range(1,4):
                            for mu in range(nu):
                                plaq += self.getLocalPlaquette(x,y,z,t,mu,nu)
        return plaq/(self.Ns**3*self.Nt*18)


    def getLinkTrace(self):
        """ Calculate <tr U>. """
        linkTrace = 0.
        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(4):
                            linkTrace += self.getLink(x,y,z,t,mu).trace().real
        linkTrace /= (self.Ns**3*self.Nt*4*3)
        return linkTrace
