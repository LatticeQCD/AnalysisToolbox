# 
# main_nuclearMatter.py                                                               
# 
# D. Clarke
# 
# This example is just to showcase some of the useful features of
# the AnalysisToolbox to facilitate doing back-of-the-envelope calculations 
# in a nuclear physics context
#

import numpy as np
from latqcdtools.physics.constants import M_proton_phys, M_neutron_phys
import latqcdtools.base.logger as logger

#
# Here we estimate some useful scales when trying to get an intuition
# about nuclear matter density
#
M_n  = ( M_proton_phys(units="GeV")+M_neutron_phys(units="GeV") )/2 # nucleon mass
r_n  = 0.84 # approximate nucleon radius (Recent results seem to favor 0.84 over 0.87)
V_n  = (4/3)*np.pi*r_n**3
pack = np.pi/3/np.sqrt(2) # close packing of equal spheres (Gauss)
n_0  = 0.15 # nuclear saturation density inferred from experiment (https://doi.org/10.1103/PhysRevC.102.044321)

eps = M_n/V_n

logger.info('Energy density at nuclear saturation: ',round(n_0*M_n,2) ,'[GeV/fm^3]')
logger.info('Energy density at close-packing limit:',round(pack*eps,2),'[GeV/fm^3]')
logger.info('Energy density within a nucleon:      ',round(eps,2)     ,'[GeV/fm^3]')

