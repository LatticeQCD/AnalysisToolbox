# 
# testCheck.py                                                               
# 
# Make sure the internal consistency methods do what they are supposed to. 
# 

import numpy as np
from latqcdtools.base.check import checkType, checkEqualLengths
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testCheck():

    checkType(int,test=1) 
    checkType(bool,test=True) 
    checkType(list,test=[1,'two',3])
    checkType("array",test=[1,'two',3])
    checkType("array",test=np.array([1,'two',3]))
    checkType("array",test=(None,None,1))
    checkType("scalar",test=1)
    checkType("scalar",test=1.)
    checkType("scalar",test=3+1j)
    checkType("real",test=3)
    checkType("real",test=np.complex128(3))

    checkEqualLengths([1,1,1],[2,2,2],[3,3,3],[4,4,4])

    logger.TBPass('All tests passed.')


if __name__ == '__main__':
    testCheck()
