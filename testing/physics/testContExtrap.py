# 
# testContExtrap.py                                                               
# 
# D. Clarke
# 
# Test of methods carrying out continuum-limit extrapolations.
# 

from latqcdtools.base.check import print_results
from latqcdtools.physics.continuumExtrap import extrapolate_from_a

PREC = 1e-7

a        = [0.09, 0.12, 0.15]
a_mu     = [-3.83725749e-01, -2.50780435e-01, -1.51850559e-01]
a_mu_err = [7.05719861e-03, 1.60531523e-02, 9.46881142e-03]

result, result_err, chidof = extrapolate_from_a(a,a_mu,a_mu_err,show_results=True,plot_results=False,obsName='$a_\mu^{\\rm disc.}$')

REFresult     = [-0.51180259, 16.18682616]
REFresult_err = [0.02250462, 1.51298996]
REFchidof     = 3.4141971374256856

print_results(result,REFresult,result_err,REFresult_err,text='result',prec=PREC)
print_results(chidof,REFchidof,text='chi^2/d.o.f.',prec=PREC)
