#!/usr/bin/env python3

import argparse
import latqcdtools.experimental.solve as slv
import latqcdtools.experimental.tools as tls
import latqcdtools.experimental.scales_hisq as sq


def get_T_GeV(beta, nt, year, suppress_warnings=False):
    return 1. / (nt * sq.a_r1_invGeV(beta, year, suppress_warnings))


def search_beta(T_GeV, nt, year):
    return slv.solve(get_T_GeV, T_GeV, 1, 15, 10e-12, nt, year, True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--T", type=float)
    parser.add_argument("--nt", type=int)
    parser.add_argument("--beta", type=float)
    parser.add_argument("--year", type=str, help="use the fit parameters from this year", required=True)
    args = parser.parse_args()

    if args.beta and not args.nt:
        beta = args.beta
        print("Beta: %f a: %f 1/GeV a^-1: %f GeV" % (beta, sq.a_r1_invGeV(beta, args.year), 1 / sq.a_r1_invGeV(beta, args.year)))

    elif args.T and args.nt:
        beta = search_beta(args.T, args.nt, args.year)
        print("Beta: %f \t T: %fGeV \t nt: %i \t a: %f 1/GeV \t Tc: 0.154 GeV"
              % (beta, args.T, args.nt, sq.a_r1_invGeV(beta, args.year)))

    elif args.beta and args.nt:
        print("Beta: %f \t T: %fGeV \t nt: %i \t a: %f 1/GeV \t Tc: 0.154 GeV"
              % (args.beta, get_T_GeV(args.beta, args.nt, args.year), args.nt, sq.a_r1_invGeV(args.beta, args.year)))
    else:
        print("Usage " + tls.sys.argv[0] + " --nt nt --beta beta --T T")
