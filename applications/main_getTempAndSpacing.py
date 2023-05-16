#!/bin/python3

#
# main_getTempAndSpacing.py
#
# D. Clarke
#
# Given an Nt, beta, and a reference scale, get a and T in physical units.
#

import argparse
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.base.utilities import getArgs
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

parser = argparse.ArgumentParser(description='Compute a and T.')

parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
parser.add_argument('--Ns', dest='Ns', type=int, default=None, help='spatial extension')
parser.add_argument('--beta', dest='beta', type=float, required=True, help='bare coupling const')
parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')
parser.add_argument('--year', default=2021, help='select the fit parameters by specifying the year of the paper')

args = getArgs(parser)

if args.Ns is None:
    Ns = 3*args.Nt
else:
    Ns = args.Ns

lp = latticeParams(Ns, args.Nt, args.beta, mass1=None, mass2=None, scaleType=args.scale, paramYear=args.year)

lp.paramSummary()
