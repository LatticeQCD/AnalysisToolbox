# 
# SU3.py                                                               
#
# D. Clarke
#
# Implementation of SU3 gauge group in Python. We try to implement various functions as static methods 
# so that we can compile them with numba to speed them up. We also try to leverage numpy.
# 

import numpy as np
from latqcdtools.math.math import SUN
from latqcdtools.base.speedify import compile


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


class SU3(SUN):

    def __new__(cls, mat=None):
        return super(SU3, cls).__new__(cls, N=3, mat=mat)

    def __repr__(self) -> str:
        return "SU(3)"

    def su3unitarize(self):
        """ 
        Project to SU(3) using information from the first two rows. 
        """
        fastUnitarize(self)

    def setToRandom(self):
        """ 
        Turn into a randomly chosen SU(3) matrix. 
        """
        fastRandomize(self)
        self.su3unitarize()
