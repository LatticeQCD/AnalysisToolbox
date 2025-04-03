# 
# SU3.py                                                               
#
# D. Clarke
#
# Implementation of SU3 gauge group in Python. We try to implement various functions as static methods 
# so that we can compile them with numba to speed them up. We also try to leverage numpy.
# 

import numpy as np
from numpy.linalg import det
import latqcdtools.base.logger as logger
from latqcdtools.math.math import id_3, ze_3, isUnitary, isSpecial, dagger
from latqcdtools.base.speedify import compile
from latqcdtools.base.check import checkType


# Eventually we would like to use default_rng here too, but it doesn't compile straightforwardly
# using numba. Will have to think about it.
rng = np.random


@compile
def fastUnitarize(self):

    # Normalize row 0.
    norm = np.sqrt( abs(self[0,0])**2 + abs(self[0,1])**2 + abs(self[0,2])**2 )
    self[0,0] /= norm
    self[0,1] /= norm
    self[0,2] /= norm

    # Get row 1 from ON projection on row 0.
    Cre =     self[1,0].real*self[0,0].real + self[1,0].imag*self[0,0].imag \
            + self[1,1].real*self[0,1].real + self[1,1].imag*self[0,1].imag \
            + self[1,2].real*self[0,2].real + self[1,2].imag*self[0,2].imag
    Cim =     self[1,0].imag*self[0,0].real - self[1,0].real*self[0,0].imag \
            + self[1,1].imag*self[0,1].real - self[1,1].real*self[0,1].imag \
            + self[1,2].imag*self[0,2].real - self[1,2].real*self[0,2].imag
    self[1,0] -= complex( Cre*self[0,0].real - Cim*self[0,0].imag, Cre*self[0,0].imag + Cim*self[0,0].real )
    self[1,1] -= complex( Cre*self[0,1].real - Cim*self[0,1].imag, Cre*self[0,1].imag + Cim*self[0,1].real )
    self[1,2] -= complex( Cre*self[0,2].real - Cim*self[0,2].imag, Cre*self[0,2].imag + Cim*self[0,2].real )

    # Normalize row 1.
    norm = np.sqrt( abs(self[1,0])**2 + abs(self[1,1])**2 + abs(self[1,2])**2 )
    self[1,0] /= norm
    self[1,1] /= norm
    self[1,2] /= norm

    # Row 2 is vector product of rows 0 and 1.
    self[2,0] = complex( ( self[0,1].real*self[1,2].real - self[0,1].imag*self[1,2].imag ) - ( self[0,2].real*self[1,1].real - self[0,2].imag*self[1,1].imag ) ,
                        -( self[0,1].imag*self[1,2].real + self[0,1].real*self[1,2].imag ) + ( self[0,2].imag*self[1,1].real + self[0,2].real*self[1,1].imag ) )

    self[2,1] = complex( ( self[0,2].real*self[1,0].real - self[0,2].imag*self[1,0].imag ) - ( self[0,0].real*self[1,2].real - self[0,0].imag*self[1,2].imag ),
                        -( self[0,2].imag*self[1,0].real + self[0,2].real*self[1,0].imag ) + ( self[0,0].imag*self[1,2].real + self[0,0].real*self[1,2].imag ) )

    self[2,2] = complex( ( self[0,0].real*self[1,1].real - self[0,0].imag*self[1,1].imag ) - ( self[0,1].real*self[1,0].real - self[0,1].imag*self[1,0].imag ),
                        -( self[0,0].imag*self[1,1].real + self[0,0].real*self[1,1].imag ) + ( self[0,1].imag*self[1,0].real + self[0,1].real*self[1,0].imag ) )


@compile
def fastRandomize(self):
    for i in range(3):
        for j in range(3):
            self[i,j] = complex( 1 - 2*rng.uniform(0,1), 1 - 2*rng.uniform(0,1) )


class SU3(np.ndarray):

    """ 
    A member of the Lie group SU(3). Implemented as a subclass of the np.ndarray class. This gives us access already
    to all the nice features of np.ndarray and lets us leverage the speed of numpy.
        g.trace()
        g.det()
        g.dagger()
        g[i,j], which can be used to access and assign
        g + h
        g*h = g@h
        2*g
    """

    def __new__(cls, mat=None):
        if mat is None:
            mat = ze_3 
        obj = np.asarray(mat, dtype=complex)
        if obj.shape != (3, 3):
            logger.TBRaise("SU(3) matrices must have shape (3,3).")
        return np.copy(obj).view(cls)


    def __repr__(self) -> str:
        return "SU(3)"


    def __mul__(self, other):
        # Perform matrix multiplication instead of element-wise multiplication
        return np.dot(self, other)


    def __pow__(self,power):
        # Perform matrix power instead of element-wise power
        checkType(int,power=power)
        return np.linalg.matrix_power(self, power)


    def trace(self):
        """ 
        Trace. In np.matrix, this returns a 2d object for some reason. 
        """
        return super().trace().item()


    def dagger(self):
        return dagger(self) 


    def det(self):
        """ 
        Determinant. 
        """
        return det(self)


    def isSU3(self) -> bool:
        """ 
        Check that I have det=1 and am unitary. 
        """
        if not (isSpecial(self) and isUnitary(self)):
            return False
        return True


    def su3unitarize(self):
        """ 
        Project to SU(3) using information from the first two rows. 
        """
        fastUnitarize(self)


    def setToMatrix(self,other):
        """ 
        Turn into RHS link. 
        """
        np.copyto(self,np.asarray(other,dtype=complex))


    def setToZero(self):
        """ 
        Turn into zero matrix. 
        """
        self.setToMatrix(ze_3)


    def setToIdentity(self):
        """ 
        Turn into identity matrix. 
        """
        self.setToMatrix(id_3) 


    def setToRandom(self):
        """ 
        Turn into a randomly chosen SU(3) matrix. 
        """
        fastRandomize(self)
        self.su3unitarize()
