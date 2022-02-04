#!/usr/bin/env python3

import numpy as np
import argparse, sys
from latqcdtools.solve import *
from latqcdtools.scales_quenched import *

parser = argparse.ArgumentParser()

parser.add_argument("--T", type=float) 
parser.add_argument("--nt", type=int) 
parser.add_argument("--beta", type=float) 
args = parser.parse_args()


def get_T_GeV(beta, nt):
    return 1./(nt*a_r0_invGeV(beta))

def get_Tc_GeV():
    r0Tc = 0.7457
    return r0Tc / r0_phys_GeV


def search_beta(T_GeV, nt):
    return solve(get_T_GeV, T_GeV, 5.7, 7.8, 1e-12, nt)


Tc=get_Tc_GeV()
nt=args.nt

if args.T and args.nt:
    T=args.T
    beta = search_beta(args.T, args.nt)
elif args.beta and args.nt:
    T=get_T_GeV(args.beta, args.nt)
    beta = args.beta
else:
    print("Usage " + sys.argv[0] + " --nt nt --beta beta --T T")
    sys.exit(-1)

print("Beta: %f \t T: %fGeV\t T/T_c: %f \t nt: %i \tr_0/a: %f\t a: %f 1/GeV\ta: %f fm\t ""1/a: %f GeV"
        "\t T_c: %fGeV" % (beta, T, T/Tc, nt, np.exp(ln_r0(beta)), a_r0_invGeV(beta), a_r0_fm(beta),
            1/a_r0_invGeV(beta), Tc))