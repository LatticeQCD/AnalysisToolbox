# 
# testStatPhys.py                                                               
# 
# D. Clarke
# 
# Test some of the statistical physics methods.
#

from latqcdtools.physics.statisticalPhysics import O2_3d, O3_3d, O4_3d, Z2_4d, Z2_3d, Z2_2d
import latqcdtools.base.logger as logger

univ = O2_3d()
univ.hyperscalingCheck()
univ = O3_3d()
univ.hyperscalingCheck(tol=1e-2)
univ = O4_3d()
univ.hyperscalingCheck()
univ = Z2_4d()
univ.hyperscalingCheck()
univ = Z2_3d()
univ.hyperscalingCheck()
univ = Z2_2d()
univ.hyperscalingCheck()

logger.TBPass("All hyperscaling checks passed.")
