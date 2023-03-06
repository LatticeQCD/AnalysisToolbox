#!/bin/python3

#
# main_getBeta.py
#
# D. Clarke
#
# Given an Nt and T, get beta.
#

import argparse
from latqcdtools.physics.referenceScales import r0_div_a,r0_hQCD_2014
from latqcdtools.base.utilities import getArgs
from latqcdtools.math.optimize import persistentSolve
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

parser = argparse.ArgumentParser(description='Compute a and T.')

parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')
parser.add_argument('--T', dest='T', required=True, type=int, help='temperature in [MeV]')

args = getArgs(parser)

scale      = args.scale
T          = args.T
beta_guess = 6

# Implement target_T - T as the LHS of the solver.
if scale == 'r0':
    def LHS(beta):
        return r0_div_a(beta)/r0_hQCD_2014("MeVinv")/args.Nt - T
else:
    logger.TBError('Beta extraction not yet implemented for scale',scale)

print('   scale =',scale)
print('       T = ',T,'[MeV]')
print('      Nt = ',args.Nt)

beta_solution = persistentSolve(LHS,beta_guess)

print('beta_sol =',round(float(beta_solution),3))