# 
# testSU3.py                                                               
# 
# D. Clarke
# 
# Test some of the basic methods for SU3 matrices.
# 

from latqcdtools.math.SU3 import SU3, id_3
import latqcdtools.base.logger as logger

g = SU3()
h = SU3()

ltest = True

# Assignment and comparison tests
g[0,0]=2
h[0,0]=2

if not g.isEqualTo(h):
   ltest = False
   logger.TBFail('Assignment and comparison.')

# Trace
if g.trace() != complex(4):
    ltest = False
    logger.TBFail('Trace.')

# Making a random matrix
g.setToRandom()
if not g.isSU3():
    ltest = False
    logger.TBFail('Set to random.')

g.setToIdentity()
if not g.isEqualTo(id_3):
    ltest = False
    logger.TBFail('Set to identity.')

if ltest:
    logger.TBPass('SU3 tests passed.')