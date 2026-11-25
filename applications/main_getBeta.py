#!/usr/bin/env python3

#
# main_getBeta.py
#
# D. Clarke, K. Ebira, Google Antigravity
#
# Given an Nt and T, get beta.
#

import argparse
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.physics.referenceScales import r0_div_a, a_times_fk, a_div_r1, ignoreBetaRange
from latqcdtools.physics.constants import convert
from latqcdtools.base.utilities import getArgs, unvector
from latqcdtools.math.optimize import persistentSolve
import latqcdtools.base.logger as logger


def main():
    logger.set_log_level('INFO')
    ignoreBetaRange()

    parser = argparse.ArgumentParser(description='Given Nt and T, solve for bare coupling beta.')
    parser.add_argument('--Nt', dest='Nt', type=int, default=8, help='euclidean time extension')
    parser.add_argument('--Nf', dest='Nf', type=str, default=None,
                        help='number flavors, e.g. 21 or None (defaults to None for r0, 21 otherwise)')
    parser.add_argument('--scale', dest='scale', required=True, help='reference scale (r0, r1, or fk)')
    parser.add_argument('--T', dest='T', required=True, type=float, help='temperature in [MeV]')
    parser.add_argument('--paramYear', dest='paramYear', type=int, default=None, help='year for a(beta) parameterization')
    parser.add_argument('--scaleYear', dest='scaleYear', type=int, default=None, help='year for scale in physical units')
    parser.add_argument('--physWorld', dest='physWorld', type=str, default=None,
                        help="world for physical scale in fm/MeV (for fk: 'nature', 'Nf2', 'Nf21', or 'Nf211')")

    args = getArgs(parser)

    scale      = args.scale
    T          = args.T
    Nt         = args.Nt
    Nf         = args.Nf
    scaleYear  = args.scaleYear
    paramYear  = args.paramYear
    physWorld  = args.physWorld
    beta_guess = 6.0

    if scale not in ('r0', 'r1', 'fk'):
        logger.TBRaise('Beta extraction not yet implemented for scale', scale)

    if Nf is None:
        Nf = None if scale == 'r0' else '21'
    elif Nf.strip().lower() in ('none', '0'):
        Nf = None

    # Reuse latticeParams as the single source of truth for Nf validation, world resolution, and physical scales.
    lp = latticeParams(None, Nt, beta_guess, Nf=Nf, scaleType=scale,
                       paramYear=paramYear, scaleYear=scaleYear, physWorld=physWorld)

    # Implement T(beta) - target_T as the LHS of the solver.
    if scale == 'r0':
        r0_mevinv = convert(lp.r0, 'fm', 'MeVinv')
        def LHS(beta):
            return r0_div_a(beta, lp.year) / r0_mevinv / Nt - T
    elif scale == 'r1':
        r1_mevinv = convert(lp.r1, 'fm', 'MeVinv')
        def LHS(beta):
            return 1.0 / (a_div_r1(beta, lp.year) * r1_mevinv * Nt) - T
    else:
        fk_mev = lp.fK
        def LHS(beta):
            return fk_mev / a_times_fk(beta, lp.year) / Nt - T

    logger.info('     scale =', scale)
    logger.info('        Nf =', lp.Nf)
    logger.info('    aWorld =', lp.aWorld)
    logger.info(' physWorld =', lp.physWorld)
    logger.info(' scaleYear =', lp.scaleYear)
    logger.info(' paramYear =', lp.year)
    logger.info('         T =', T, '[MeV]')
    logger.info('        Nt =', Nt)

    beta_solution = persistentSolve(LHS, beta_guess)

    logger.info('  beta_sol =', round(unvector(beta_solution), 3))


if __name__ == '__main__':
    main()
