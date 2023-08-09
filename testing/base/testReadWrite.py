# 
# testReadWrite.py                                                               
# 
# D. Clarke
# 
# Tests for read and write methods.
#
import argparse
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import writeTable
from latqcdtools.base.utilities import getArgs

logger.set_log_level('INFO')

def testReadWrite():

    xdata = np.linspace(0,1,101)
    ydata = xdata**2

    writeTable("test.d",xdata,ydata,header=["who","bastank"])

    logger.TBPass("I didn't encounter an error.")

if __name__ == '__main__':
    testReadWrite()