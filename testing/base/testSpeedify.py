# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import numpy as np
from latqcdtools.math.math import print_results
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import parallel_function_eval, parallel_reduce, compile, numbaON, getMaxThreads

logger.set_log_level('INFO')

numbaON()

@compile
def square(x):
    return x**2

testArray = range(10)

def testSpeedify():

    temp = parallel_function_eval(square,testArray)

    sum1 = 0.
    for i in testArray:
        sum1 += square(i)
    sum2 = np.sum(temp)

    print_results(sum1,sum2,text='parallel square')

    sum3 = parallel_reduce(square,testArray)
    print_results(sum2,sum3,text='parallel square reduction')

    temp = parallel_function_eval(square,testArray,nproc=1)


if __name__ == '__main__':
    testSpeedify()