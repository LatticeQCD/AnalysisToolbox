# 
# test_Lime.py                                                               
# 
# D. Clarke
# 
# Testing the Lime interface. 
# 
from latqcdtools.interfaces.lime import limeHeader, trimNull
from latqcdtools.testing import print_results,concludeTest
import latqcdtools.base.logger as logger
import numpy as np

def testLime():

    lpass = True

    testByteString = limeHeader(1,0,123,b'test') 

    lpass *= print_results( len(testByteString), 144, text="header length")

    if trimNull(testByteString) !=  b'Eg\x89\xab':
        lpass = False
        logger.TBFail('null trim')

    concludeTest(lpass)


if __name__ == '__main__':
    testLime()
