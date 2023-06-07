# 
# gauge.py                                                               
# 
# D. Clarke
# 
# A gauge field class. We try to implement various functions as static methods so that we can compile
# them with numba to speed them up.
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.math.SU3 import SU3
from latqcdtools.base.speedify import compile, parallel_reduce
from latqcdtools.base.check import checkType


@compile
def respectBCs(x,y,z,t,Ns,Nt):
    """ Assume periodic BCs. This also seems to work for negative numbers, e.g. -1 % 8 = 7. """
    return x%Ns, y%Ns, z%Ns, t%Nt


@compile
def stepUp(x,y,z,t,mu,Ns,Nt):
    X, Y, Z, T = x, y, z, t
    if mu==0:
        X += 1
    elif mu==1:
        Y += 1
    elif mu==2:
        Z += 1
    elif mu==3:
        T += 1
    return respectBCs(X, Y, Z, T, Ns, Nt)


@compile
def stepDown(x,y,z,t,mu,Ns,Nt):
    X, Y, Z, T = x, y, z, t
    if mu==0:
        X -= 1
    elif mu==1:
        Y -= 1
    elif mu==2:
        Z -= 1
    elif mu==3:
        T -= 1
    return respectBCs(X, Y, Z, T, Ns, Nt)


@compile
def ReTrABCD(A,B,C,D):
    return np.trace( A @ B @ C @ D ).real


class gaugeField:

    """ A class for SU3 gauge fields. """

    def __init__(self, Ns, Nt, nproc=1):
        """ Initialize an SU3 gaugeField.

        Args:
            Ns (int): Spatial extension of lattice. 
            Nt (int): Euclidian time extension of lattice. 
            nproc (int, 1): Number of proccesors for parallelization. Defaults to 1.
        """
        logger.details('Initialize gaugeField object...')
        checkType(Ns,int)
        checkType(Nt,int)
        checkType(nproc,int)
        self.Ns    = Ns
        self.Nt    = Nt
        self.nproc = nproc
        self.field = np.zeros(shape=(self.Nt,self.Ns,self.Ns,self.Ns,4,3,3),dtype=complex)
        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(1,4):
                            U = SU3()
                            self.field[t,z,y,x,mu] = U


    def getLink(self,x,y,z,t,mu):
        """ Link accessor. """
        return self.field[t][z][y][x][mu]


    def setLink(self, other, x, y, z, t, mu):
        """ Link setter. """
        self.field[t][z][y][x][mu] = other


    def getLocalPlaquette(self,x,y,z,t,mu,nu):
        """ Get Re tr U_{mu,nu}^[](x,y,z,t). """
        U_mu_x     = SU3()
        U_nu_xpamu = SU3()
        U_mu_xpanu = SU3()
        U_nu_x     = SU3()
        U_mu_x.setToMatrix( self.getLink(x,y,z,t,mu) )
        X, Y, Z, T = stepUp(x,y,z,t,mu,self.Ns,self.Nt)
        U_nu_xpamu.setToMatrix( self.getLink(X,Y,Z,T,nu) )
        X, Y, Z, T = stepUp(x,y,z,t,nu,self.Ns,self.Nt)
        U_mu_xpanu.setToMatrix( self.getLink(X,Y,Z,T,mu) )
        U_nu_x.setToMatrix( self.getLink(x,y,z,t,nu) )
        return ReTrABCD(U_mu_x,U_nu_xpamu,U_mu_xpanu.dagger(),U_nu_x.dagger())


    def getPlaquette(self):
        """ Calculate <Re tr U_{mu,nu}^[]>. """
        plaq = parallel_reduce(self.plaq_contrib,range(self.Nt),self.nproc)
        return plaq/(self.Ns**3*self.Nt*18)


    def plaq_contrib(self,t):
        contrib = 0
        for z in range(self.Ns):
            for y in range(self.Ns):
                for x in range(self.Ns):
                    for nu in range(1,4):
                        for mu in range(nu):
                            contrib += self.getLocalPlaquette(x,y,z,t,mu,nu)
        return contrib


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
