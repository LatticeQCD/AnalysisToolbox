#!/usr/bin/env python3

#
# main_getTempAndSpacing.py
#
# D. Clarke, K. Ebira, Google Antigravity
#
# Given an Nt, beta, and a reference scale, get a and T in physical units.
#

import argparse
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.base.utilities import getArgs
import latqcdtools.base.logger as logger


def main():
    logger.set_log_level('INFO')

    parser = argparse.ArgumentParser(description='Compute a and T.')
    parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
    parser.add_argument('--Ns', dest='Ns', type=int, default=None, help='spatial extension')
    parser.add_argument('--Nf', dest='Nf', type=str, default=None,
                        help='number flavors (e.g. 21, 111, or None; defaults to None for r0, 21 otherwise)')
    parser.add_argument('--beta', dest='beta', type=float, required=True, help='bare coupling const')
    parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')
    parser.add_argument('--paramYear', dest='paramYear', type=int, default=None, help='year for a(beta) parameterization')
    parser.add_argument('--scaleYear', dest='scaleYear', type=int, default=None, help='year for scale in physical units')
    parser.add_argument('--physWorld', dest='physWorld', type=str, default=None,
                        help="world for physical scale in fm/MeV (for fk: 'nature', 'Nf2', 'Nf21', or 'Nf211')")

    args = getArgs(parser)

    scale      = args.scale
    Nt         = args.Nt
    Nf         = args.Nf
    beta       = args.beta
    scaleYear  = args.scaleYear
    paramYear  = args.paramYear
    physWorld  = args.physWorld

    if Nf is None:
        Nf = None if scale == 'r0' else '21'
    elif Nf.strip().lower() in ('none', '0'):
        Nf = None

    if args.Ns is None:
        Ns = 3 * args.Nt
    else:
        Ns = args.Ns

    lp = latticeParams(Ns, Nt, beta, Nf=Nf, mass1=None, mass2=None, scaleType=scale, paramYear=paramYear,
                       scaleYear=scaleYear, physWorld=physWorld)

    lp.paramSummary()


if __name__ == '__main__':
    main()
