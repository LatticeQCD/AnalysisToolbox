# 
# testCheck.py                                                               
# 
# Make sure the internal consistency methods do what they are supposed to. 
# 

import numpy as np
from latqcdtools.base.check import checkType, checkEqualLengths
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def func1(x):
    return np.log(x)

def func2(x,p):
    A, B, C = p
    return A*( 1 + B*x**2 + C*x**2*np.log(x) )

def func3(x,p,opt):
    A, B, C = p
    if opt:
        logger.info('I found my optional.')
    return A*( 1 + B*x**2 + C*x**2*np.log(x) )


def testCheck():

    checkType(1,int) 
    checkType(True,bool) 
    checkType([1,'two',3],list)
    checkType([1,'two',3],'array')
    checkType(np.array([1,'two',3]),'array')
    checkType((None,None,1),'array')

    checkEqualLengths([1,1,1],[2,2,2],[3,3,3],[4,4,4])

    logger.TBPass('All tests passed.')


if __name__ == '__main__':
    testCheck()
