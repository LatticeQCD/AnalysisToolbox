#!/bin/python3

#
# main_getBeta.py
#
# D. Clarke
#
# Given an Nt and T, get beta.
#

import argparse
from latqcdtools.physics.referenceScales import r0_div_a, r0_hQCD_2014, a_times_fk, fk_phys, ignoreBetaRange, CY_FK_PHYS, CY_A_TIMES_FK, CY_A_DIV_R0
from latqcdtools.base.utilities import getArgs
from latqcdtools.math.optimize import persistentSolve
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')
ignoreBetaRange()

parser = argparse.ArgumentParser(description='Compute a and T.')

parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')
parser.add_argument('--T', dest='T', required=True, type=int, help='temperature in [MeV]')
parser.add_argument('--paramYear', dest='paramYear', default=None, help='year for a(beta) parameterization')
parser.add_argument('--scaleYear', dest='scaleYear', default=None, help='year for scale in physical units')

args = getArgs(parser)

scale      = args.scale
T          = args.T
Nt         = args.Nt
scaleYear  = args.scaleYear
paramYear  = args.paramYear
beta_guess = 6

# Implement target_T - T as the LHS of the solver.
if scale == 'r0':
    if paramYear is None:
        paramYear = CY_A_DIV_R0
    def LHS(beta):
        return r0_div_a(beta,paramYear)/r0_hQCD_2014("MeVinv")/Nt - T
elif scale == 'fk':
    if scaleYear is None:
        scaleYear = CY_FK_PHYS
    if paramYear is None:
        paramYear = CY_A_TIMES_FK 
    def LHS(beta):
        return fk_phys(scaleYear,units="MeV")/a_times_fk(beta,paramYear)/Nt - T
else:
    logger.TBError('Beta extraction not yet implemented for scale',scale)

logger.info('     scale =',scale)
logger.info(' scaleYear =',scaleYear)
logger.info(' paramYear =',paramYear)
logger.info('         T =',T,'[MeV]')
logger.info('        Nt =',Nt)

beta_solution = persistentSolve(LHS,beta_guess)

logger.info('  beta_sol =',round(float(beta_solution),3))