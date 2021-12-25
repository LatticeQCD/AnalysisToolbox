#
# gaudifTest.py
#
# D. Clarke
#
# Simple test for the Gaussian difference test (Z-test).
#
from latqcdtools.statistics.statistics import gaudif
import latqcdtools.base.logger as logger

eps=1e-7

# Some measurements and error
x1=0.4655
e1=0.0088
x2=0.501
e2=0.045
x3=0.480
e3=0.023
x4=0.4745
e4=0.0063

q12=gaudif(x1,e1,x2,e2)
q13=gaudif(x1,e1,x3,e3)
q14=gaudif(x1,e1,x4,e4)

# Results produced by software of "Markov Chain Monte Carlo Simulations and Their Statistical Analysis, World
# Scientific, Singapore, 2004.
q12control=0.4387984
q13control=0.5559897
q14control=0.4056413

# The simple test.
lerror=False
if abs(q12-q12control) > eps:
  lerror=True
if abs(q13-q13control) > eps:
  lerror=True
if abs(q14-q14control) > eps:
  lerror=True

if lerror:
  logger.TBError("At least one test failed!")
else:
  logger.TBPass("All tests passed!")
