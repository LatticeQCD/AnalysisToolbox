# 
# testGauge.py                                                               
# 
# D. Clarke
# 
# Quick tests of the gaugeField class.
#

import latqcdtools.base.logger as logger
from latqcdtools.math.math import rel_check
from latqcdtools.math.SU3 import SU3
from latqcdtools.physics.gauge import gaugeField
from latqcdtools.testing import concludeTest


def testGauge():

    lpass = True

    gauge = gaugeField(8,4)

    g = SU3()

    g.setToRandom()
    gauge.setLink(g,0,0,0,0,0)

    # Boundary condidition tests
    h = gauge.getLink(8,0,0,0,0)
    if not rel_check(g,h):
        lpass = False
        logger.TBFail('X direction BC.')
    h = gauge.getLink(0,8,0,0,0)
    if not rel_check(g,h):
        lpass = False
        logger.TBFail('Y direction BC.')
    h = gauge.getLink(0,0,8,0,0)
    if not rel_check(g,h):
        lpass = False
        logger.TBFail('Z direction BC.')
    h = gauge.getLink(0,0,0,4,0)
    if not rel_check(g,h):
        lpass = False
        logger.TBFail('T direction BC.')

    gauge.makeCold()

    # Plaquette normalization test
    plaq = gauge.getPlaquette()
    if not rel_check(plaq,1):
        lpass = False
        logger.TBFail('Plaquette normalization. plaq =',plaq)

    concludeTest(lpass) 


if __name__ == '__main__':
    testGauge()