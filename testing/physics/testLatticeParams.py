# 
# testLatticeParams.py                                                               
# 
# D. Clarke
# 
# Some tests of the lattice parameter class using 2021 scales.
#

from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.physics.HotQCD import HotQCDParams
from latqcdtools.base.check import print_results

import params

lp = latticeParams(params.Ns, params.Nt, params.cbeta, params.cml, params.cms)
lp.paramSummary()
a = lp.geta()
T = lp.getT()
print_results(a,0.14027137436021253,text='fK test, a')
print_results(T,175.84394864955703,text='fK test, T')
del lp

lp = HotQCDParams(params.Ns, params.Nt, params.cbeta, params.cml, params.cms, Nf='3')
lp.paramSummary()
print(lp.getcparams())
del lp

import params_zeroTemp

lp = latticeParams(params_zeroTemp.Ns, params_zeroTemp.Nt, params_zeroTemp.cbeta, scaleType='r0')
lp.paramSummary()
a = lp.geta()
T = lp.getT()
print_results(a,0.04195979097233082,text='r0 test, a')
print_results(T,0.0,text='r0 test, T')
del lp
