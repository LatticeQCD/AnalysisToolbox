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

parser = argparse.ArgumentParser(description='Compute a and T.')

parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
parser.add_argument('--beta', dest='beta', type=float, default=6.285, help='bare coupling const')
parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')

args = getArgs(parser)

lp = latticeParams(None, args.Nt, args.beta, None, None, args.scale)

lp.paramSummary()
