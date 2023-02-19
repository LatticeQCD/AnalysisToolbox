# 
# testContExtrap.py                                                               
# 
# D. Clarke
# 
# Test of methods carrying out continuum-limit extrapolations.
# 

import numpy as np
from latqcdtools.base.check import print_results
from latqcdtools.physics.continuumExtrap import extrapolate_from_a
from latqcdtools.physics.constants import MeV_to_fminv

PREC = 1e-7

a        = np.array( [0.09, 0.12, 0.15] )
a_mu     = [-3.83725749e-01, -2.50780435e-01, -1.51850559e-01]
a_mu_err = [7.05719861e-03, 1.60531523e-02, 9.46881142e-03]

result, result_err, chidof = extrapolate_from_a(a,a_mu,a_mu_err,show_results=True,plot_results=False)

REFresult     = [-0.51180259, 16.18682616]
REFresult_err = [0.01217944250306863, 0.8188262377288199]
REFchidof     = 3.4141971374256856

print_results(result,REFresult,result_err,REFresult_err,text='simple O(a^2)',prec=PREC)
print_results(chidof,REFchidof,text='O(a^2) chi^2/d.o.f.',prec=PREC)

# When doing a continuum extrapolation with priors, you would like to have a rough idea of the relative sizes of the
# priors (fit coefficients). This is easier to do if the independent variable (the lattice spacing) is unitless.
# Here we shift a --> a*Lambda, where Lambda ~ Lambda_QCD ~ 500 MeV to accomplish this.

lam = MeV_to_fminv(500)
a *= lam

result, result_err, chidof = extrapolate_from_a(a,a_mu,a_mu_err,show_results=True,plot_results=False,order=2,
                                                prior=[-0.51180259,1,0.0],prior_err=[0.02250462,1,0.01])

REFresult     = [-5.10181392e-01,  2.50117866e+00,  4.10582713e-06]
REFresult_err = [0.06333996, 1.59019596, 8.10411157]
REFchidof     = 1.8992268798553928


print_results(result,REFresult,result_err,REFresult_err,text='O(a^4) with prior',prec=PREC)
print_results(chidof,REFchidof,text='O(a^4) chi^2/d.o.f.',prec=PREC)
