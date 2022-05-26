# 
# testStatPhys.py                                                               
# 
# D. Clarke
# 
# Test some of the statistical physics methods.
#

from latqcdtools.physics.statisticalPhysics import O2_3d, O4_3d, Z2_3d, Z2_2d
import latqcdtools.base.logger as logger

univ = O2_3d()
univ.exponentSummary()
univ.hyperscalingCheck()
univ = O4_3d()
univ.exponentSummary()
univ.hyperscalingCheck()
univ = Z2_3d()
univ.exponentSummary()
univ.hyperscalingCheck()
univ = Z2_2d()
univ.exponentSummary()
univ.hyperscalingCheck()

logger.TBPass("All hyperscaling checks passed.")
