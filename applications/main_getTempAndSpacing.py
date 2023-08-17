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
from latqcdtools.physics.referenceScales import CY_param, CY_phys 
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

parser = argparse.ArgumentParser(description='Compute a and T.')

parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
parser.add_argument('--Ns', dest='Ns', type=int, default=None, help='spatial extension')
parser.add_argument('--beta', dest='beta', type=float, required=True, help='bare coupling const')
parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')
parser.add_argument('--paramYear', dest='paramYear', default=None, help='year for a(beta) parameterization')
parser.add_argument('--scaleYear', dest='scaleYear', default=None, help='year for scale in physical units')

args = getArgs(parser)

scale      = args.scale
Nt         = args.Nt
beta       = args.beta
scaleYear  = args.scaleYear
paramYear  = args.paramYear

if args.Ns is None:
    Ns = 3*args.Nt
else:
    Ns = args.Ns

if paramYear is None:
    paramYear = CY_param[scale]
if scaleYear is None:
    scaleYear = CY_phys[scale]

lp = latticeParams(Ns, Nt, beta, mass1=None, mass2=None, scaleType=scale, paramYear=paramYear, scaleYear=scaleYear)

lp.paramSummary()
