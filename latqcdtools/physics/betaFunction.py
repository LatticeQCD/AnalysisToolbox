# 
# betaFunction.py                                                               
# 
# D. Clarke
# 
# Python implementation of the QCD beta function and its coefficients.
# 
import numpy as np


def b0(nf):
    return (11 - 2*nf/3) / (4*np.pi)**2


def b1(nf):
    return (102 - 38*nf/3) / (4*np.pi)**4


def beta_func(beta,nf=3):
    """ Asymptotic scaling relation up to two loops. Default nf=3. """
    return (b0(nf)*10/beta)**( -b1(nf)/(2*b0(nf)**2) ) * np.exp( -beta/(20*b0(nf)) )