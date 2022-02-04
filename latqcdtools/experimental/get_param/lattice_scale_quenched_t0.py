#!/usr/bin/env python3

import argparse, sys
from latqcdtools.solve import *
from latqcdtools.scales_quenched import *


parser = argparse.ArgumentParser()

parser.add_argument("--T", type=float) 
parser.add_argument("--nt", type=int) 
parser.add_argument("--beta", type=float) 
args = parser.parse_args()


def get_T_GeV(beta, nt):
    return 1./(nt*a_t0_invGeV(beta))

def get_Tc_GeV():
    return 0.2489 / sqrtt0_phys

def search_beta(T_GeV, nt):
    return solve(get_T_GeV, T_GeV, 5.7, 8.0, 10e-12, nt)


if args.T and args.nt:
    beta = search_beta(args.T, args.nt)
    print("Beta: %f \t T: %fGeV \t T/T_c: %f \t nt: %i \t a: %f 1/GeV \t T_c: %fGeV"
            % (beta, args.T, args.T/get_Tc_GeV(), args.nt, a_t0_invGeV(beta), get_Tc_GeV()))

elif args.beta and args.nt:
    print("Beta: %f \t T: %fGeV \t T/T_c: %f \t nt: %i \t a: %f 1/GeV \t T_c: %fGeV"
            % (args.beta, get_T_GeV(args.beta, args.nt),
                get_T_GeV(args.beta, args.nt)/get_Tc_GeV(), args.nt,
                a_t0_invGeV(args.beta), get_Tc_GeV()))
else:
    print("Usage " + sys.argv[0] + "--nt nt --beta beta --T T")