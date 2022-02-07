#!/usr/bin/env python3

import argparse
import latqcdtools.experimental.solve as solve
import latqcdtools.experimental.tools as tools
import latqcdtools.experimental.scales_hisq as sh


def get_T_GeV(beta, nt):
    return 1. / (nt * sh.a_fk_invGeV_2021(beta))


def search_beta(T_GeV, nt):
    return solve.solve(get_T_GeV, T_GeV, 5.9, 7.8, 10e-12, nt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--T", type=float)
    parser.add_argument("--nt", type=int)
    parser.add_argument("--beta", type=float)
    args = parser.parse_args()

    print("using data from 2021")

    if args.beta and not args.nt:
        beta = args.beta
        print("Beta: %f a: %f fm \t a: %f 1/GeV, a^-1: %f GeV, af_k: %f" % (beta, sh.a_fk_fm_2021(beta), sh.a_fk_invGeV_2021(beta),
                                                                            1 / sh.a_fk_invGeV_2021(beta), sh.a_times_fk_2021(beta)))

    elif args.T and args.nt:
        beta = search_beta(args.T, args.nt)
        print("Beta: %f \t T: %f GeV \t nt: %i \t a: %f fm \t a: %f 1/GeV \t Tc: 0.1565 GeV"
              % (beta, args.T, args.nt, sh.a_fk_fm_2021(beta), sh.a_fk_invGeV_2021(beta)))

    elif args.beta and args.nt:
        print("Beta: %f \t T: %f GeV \t nt: %i \t a: %f fm \t a: %f 1/GeV \t Tc: 0.1565 GeV"
              % (args.beta, get_T_GeV(args.beta, args.nt), args.nt, sh.a_fk_fm_2021(args.beta), sh.a_fk_invGeV_2021(args.beta)))
    else:
        print("Usage " + tools.sys.argv[0] + " --nt nt --beta beta --T T")