#!/usr/bin/env python3

import argparse
from latqcdtools.scales_hisq import *


parser = argparse.ArgumentParser()

parser.add_argument("beta", type=float) 
args = parser.parse_args()


def pol(beta):
    a, b, c, d = -2.3505170361260177, 47.42881022059038, -319.73724305951777, 720.7116076653417
    return a*beta**3 + b*beta**2 + c*beta + d

def naik_eps(mc):
    mc2=mc*mc;
    mc4=mc2*mc2;
    mc6=mc4*mc2;
    mc8=mc4*mc4;

    eps=(-27.0/40.0)*mc2 + (327.0/1120.0)*mc4 - (15607.0/268800.0)*mc6 - (73697.0/3942400.0)*mc8
    return eps


def amc_n2_dn1(beta):
    c0, c2, d2 = 55.80956116442332, 773754.1484925086, 5690.682250893555
    return (c0*beta_func(beta)+c2*10./beta*beta_func(beta)**3) / \
            (1+d2*10./beta*beta_func(beta)**2)

def amc(beta):
    if beta < 6.664:
        return pol(beta)
    else:
        return amc_n2_dn1(beta)

amc = amc(args.beta)
print("Beta: %f \t am_c: %f naik epsilon: %f" % (args.beta, amc, naik_eps(amc)))



    
    

