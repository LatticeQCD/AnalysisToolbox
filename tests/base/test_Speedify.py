# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import numpy as np
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import parallel_function_eval, parallel_reduce, compile, numbaON, \
    DEFAULTTHREADS

logger.set_log_level('INFO')

numbaON()

@compile
def square(x):
    return x**2

def power(x,i):
    return x**i

testArray = np.array(range(10))

def testSpeedify():

    lpass = True

    temp = parallel_function_eval(square,testArray,nproc=DEFAULTTHREADS)

    sum1 = 0.
    for i in testArray:
        sum1 += square(i)
    sum2 = np.sum(temp)

    lpass *= print_results(sum1,sum2,text='parallel square')

    sum3 = parallel_reduce(square,testArray,nproc=DEFAULTTHREADS)
    lpass *= print_results(sum2,sum3,text='parallel square reduction')

    temp = parallel_function_eval(power,testArray,args=(3,),nproc=DEFAULTTHREADS)

    sum1 = 0.
    for i in testArray:
        sum1 += power(i,3)
    sum2 = np.sum(temp)

    lpass *= print_results(sum1,sum2,text='parallel power')

    concludeTest(lpass)


if __name__ == '__main__':
    testSpeedify()