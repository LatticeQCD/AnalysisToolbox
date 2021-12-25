from latqcdtools.lattice_params import *

import params

lp = latticeParams(params.Ns, params.Nt, params.cbeta, params.cml, params.cms)
lp.paramSummary()
print("\n"+lp.getcparams())
del lp

import params_zeroTemp

lp = latticeParams(params_zeroTemp.Ns, params_zeroTemp.Nt, params_zeroTemp.cbeta) 
lp.paramSummary()
del lp
