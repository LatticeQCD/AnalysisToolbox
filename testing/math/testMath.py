# 
# testMath.py                                                               
# 
# D. Clarke
# 
# For testing some of the basic math functions.
#
import math
import latqcdtools.base.logger as logger
from latqcdtools.math.math import fallFactorial 
from latqcdtools.testing import print_results, concludeTest


logger.set_log_level('INFO')


def testMath():

    lpass = True

    lpass *= print_results(fallFactorial(23,23),1.*math.factorial(23),text='fallFactorial 23!')
    lpass *= print_results(fallFactorial(6,3),6*5*4,text='fallFactorial 6_fall_3')

    concludeTest(lpass)


if __name__ == '__main__':
    testMath()