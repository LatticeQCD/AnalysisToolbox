#!/usr/bin/env python3


import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("mu", type=float) 
parser.add_argument("nf", type=int) 
parser.add_argument("r0", type=float) 
args = parser.parse_args()


# four-loop running gauge coupling taken from KG. Chetyrkin et al. /Nuclear Physics B 510 (1998) 61-87
def running_coupling(mu,nf,r0):
    C3=1.202057

    beta0 = (1/4) * (11-(2/3)*nf)
    beta1 = (1/16) * (102 - (38/3) * nf)
    beta2 = (1/64) * ((2857/2) - (5033/18)*nf + (325/54)*nf**2 )
    beta3 = (1/256) * ((149753/6) + 3564*C3 + (-(1078361/162) - (6508/27)*C3)*nf + ((50065/162) + (6472/81)*C3)*nf**2 + (1093/729)*nf**3)

    b1 = beta1/beta0
    b2 = beta2/beta0
    b3 = beta3/beta0
    
    # Lambda MS taken from S. Capitani et al. [ALPHA Collaboration], Nucl. Phys. B544, 669 (1999) [hep-lat/9810063]
    lambd = 0.602/r0

    L = np.log(mu**2/lambd**2)
    lnL = np.log(L)

    b0L = 1/(beta0*L)

    a_s = (b0L - (b1*lnL)/((beta0*L)**2) + b0L**3 * (b1**2 * (lnL**2 - lnL - 1) + b2) + b0L**4 * (b1**3*(-lnL**3 + (5/2)*lnL**2 + 2*lnL - (1/2)) - 3*b1*b2*lnL + (b3/2)))*np.pi

    return a_s


print(running_coupling(args.mu,args.nf,args.r0))
