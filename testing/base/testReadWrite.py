# 
# testReadWrite.py                                                               
# 
# D. Clarke
# 
# Tests for read and write methods.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable, writeTable
from latqcdtools.testing import print_results, concludeTest 

logger.set_log_level('INFO')

def testReadWrite():

    lpass = True

    xdata = np.linspace(0,1,101)
    ydata = xdata**2

    writeTable("test.d",xdata,ydata,header=["who","bastank"])
    xtest,ytest = readTable("test.d")

    lpass *= print_results(xdata,xtest,text='x data')
    lpass *= print_results(ydata,ytest,text='y data')

    concludeTest(lpass)

if __name__ == '__main__':
    testReadWrite()