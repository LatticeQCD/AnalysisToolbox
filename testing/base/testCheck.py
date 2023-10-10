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
