# 
# testContExtrap.py                                                               
# 
# D. Clarke
# 
# Test of methods carrying out continuum-limit extrapolations.
# 

import numpy as np
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.physics.continuumExtrap import continuumExtrapolate
from latqcdtools.physics.constants import MeV_to_fminv
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testContExtrap():

    PREC = 1e-3

    lpass = True

    a         = np.array( [0.09, 0.12, 0.15] )
    a_mu      = [-3.83725749e-01, -2.50780435e-01, -1.51850559e-01]
    a_mu_err  = [7.05719861e-03, 1.60531523e-02, 9.46881142e-03]
    prior     = [-0.51180259,2,0.0]
    prior_err = [0.02250462,1,1]

    result, result_err, chidof = continuumExtrapolate(a,a_mu,a_mu_err,show_results=True)

    REFresult     = [-0.51180259, 16.18682616]
    REFresult_err = [0.01217944250306863, 0.8188262377288199]
    REFchidof     = 3.4141971374256856

    lpass *= print_results(result,REFresult,result_err,REFresult_err,text='simple O(a^2)',prec=PREC)
    lpass *= print_results(chidof,REFchidof,text='O(a^2) chi^2/d.o.f.',prec=PREC)

    # When doing a continuum extrapolation with priors, you would like to have a rough idea of the relative sizes of the
    # priors (fit coefficients). This is easier to do if the independent variable (the lattice spacing) is unitless.
    # Here we shift a --> a*Lambda, where Lambda ~ Lambda_QCD ~ 500 MeV to accomplish this.

    lam = MeV_to_fminv(500)
    a *= lam

    result, result_err, chidof, stats = continuumExtrapolate(a,a_mu,a_mu_err,show_results=True,order=2,
                                                             prior=prior,prior_err=prior_err,error_strat='hessian')

    REFresult     = [-0.51193317,  2.53542709, -0.11948167]
    REFresult_err = [0.00930215, 0.1100613,  0.10344893]
    REFchidof     = 1.2219934235330452
    REFlogGBF     = 5.36146864464183

    lpass *= print_results(result,REFresult,result_err,REFresult_err,text='O(a^4) with prior',prec=PREC)
    lpass *= print_results(chidof,REFchidof,text='O(a^4) chi^2/d.o.f.',prec=PREC)
    lpass *= print_results(stats['logGBF'],REFlogGBF,text='O(a^4) logGBF',prec=PREC)

    concludeTest(lpass)


if __name__ == '__main__':
    testContExtrap()