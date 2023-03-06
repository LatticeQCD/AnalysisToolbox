# 
# testMath.py                                                               
# 
# D. Clarke
# 
# For testing some of the basic math functions.
#
import math
import numpy as np
from latqcdtools.math.math import fallFactorial, underflowPower, underflowExp, underflowMultiply
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testMath():

    print_results(fallFactorial(23,23),1.*math.factorial(23),text='fallFactorial 23!')
    print_results(fallFactorial(6,3),6*5*4,text='fallFactorial 6_fall_3')

    smallNormal = 2**(-1022)     # Smallest 64-bit normal number

    a = np.array([1,1,smallNormal,1,1])
    b = np.array([0,0,-1000,0,0])
    c = np.array([1,1,0,1,1])

    print_results( underflowPower(a,2)     , c, text='underflowPower')
    print_results( underflowExp(b)         , c, text='underflowExp')
    print_results( underflowMultiply(a,a,a), c, text='underflowMultiply')


if __name__ == '__main__':
    testMath()