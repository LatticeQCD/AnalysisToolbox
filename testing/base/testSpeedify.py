# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import numpy as np
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import parallel_function_eval

logger.set_log_level('INFO')


def square(x):
    return x**2

testArray = range(10)

def testSpeedify():

    temp = parallel_function_eval(square,testArray,8)

    sum1 = 0.
    for i in testArray:
        sum1 += square(i)
    sum2 = np.sum(temp)

    print_results(sum1,sum2,text='parallel square test')

    temp = parallel_function_eval(square,testArray,1)


if __name__ == '__main__':
    testSpeedify()