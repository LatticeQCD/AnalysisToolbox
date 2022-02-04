#!/usr/bin/env python3

import argparse
from latqcdtools.solve import *
from latqcdtools.tools import *
from latqcdtools.scales_hisq import *


parser = argparse.ArgumentParser()

parser.add_argument("--T", type=float) 
parser.add_argument("--nt", type=int) 
parser.add_argument("--beta", type=float) 
args = parser.parse_args()

   
def get_T_GeV(beta, nt):
    return 1./(nt*a_fk_invGeV_2014(beta))

def search_beta(T_GeV, nt):
    return solve(get_T_GeV, T_GeV, 5.9, 7.8, 10e-12, nt)


if args.beta and not args.nt:
    beta = args.beta
    print("Beta: %f a: %f fm \t a: %f 1/GeV, a^-1: %f GeV, af_k: %f" % (beta, a_fk_fm_2014(beta), a_fk_invGeV_2014(beta),
        1/a_fk_invGeV_2014(beta), a_times_fk_2014(beta)))

elif args.T and args.nt:
    beta = search_beta(args.T, args.nt)
    print("Beta: %f \t T: %f GeV \t nt: %i \t a: %f fm \t a: %f 1/GeV \t Tc: 0.1565 GeV"
            % (beta, args.T, args.nt, a_fk_fm_2014(beta), a_fk_invGeV_2014(beta)))

elif args.beta and args.nt:
    print("Beta: %f \t T: %f GeV \t nt: %i \t a: %f fm \t a: %f 1/GeV \t Tc: 0.1565 GeV"
            % (args.beta, get_T_GeV(args.beta, args.nt), args.nt, a_fk_fm_2014(args.beta), a_fk_invGeV_2014(args.beta)))
else:
    print("Usage " + sys.argv[0] + " --nt nt --beta beta --T T")