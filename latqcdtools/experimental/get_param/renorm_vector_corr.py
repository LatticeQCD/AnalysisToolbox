#!/usr/bin/env python3
from latqcdtools.spline_interpolate import *
import sys
beta = float(sys.argv[1])
kappa = float(sys.argv[2])


def kappa_c(beta, betas=None, kappa_cs=None):
    if betas is None: #from Lüscher et al.
        betas = [6.0,
                6.2, 
                6.4, 
                6.8, 
                7.4, 
                8.0, 
                9.6]

        kappa_cs = [
                0.135196,
                0.135795,
                0.135720,
                0.135097,
                0.134071,
                0.133173,
                0.131448
                ]
    int_kappa_c = spline(beta, betas, kappa_cs)
    return int_kappa_c


def mq(kappa, beta, kappac=None):
    if kappac is None:
        return 1.0/(2.0*kappa)-1/(2.0*kappa_c(beta))
    else:
        return 1.0/(2.0*kappa)-1/(2.0*kappac)

g2=6.0/beta

zv=(1.-0.7663*g2+0.0488*g2**2)/(1.0-0.6369*g2)

bv=(1.0-0.6518*g2-0.1226*g2**2)/(1.0-0.8467*g2)

#print(zv, zv*(1+bv*mq(kappa,beta,0.134081)))
print(zv, zv*(1+bv*mq(kappa,beta)))
#print(zv, zv*(1+bv*10.10445/11.2))


