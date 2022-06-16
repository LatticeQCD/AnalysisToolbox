# 
# test_dataCleaner.py
# 
# D. Clarke
# 
# Some methods to test the cleanData.py module. TODO: Add test for intersectAtCol
# 


import numpy as np
from latqcdtools.base.cleanData import clipRange
import latqcdtools.base.logger as logger


testArray = np.linspace(-10,10,1000)


clippedArray = clipRange(testArray,minVal=-1,maxVal=1)
for i in range(len(clippedArray)):
    if clippedArray[i]<-1:
        logger.TBError("Found less than -1.")
    elif clippedArray[i]>1:
        logger.TBError("Found greater than 1.")


clippedArray = clipRange(testArray,minVal=-1)
for i in range(len(clippedArray)):
    if clippedArray[i]<-1:
        logger.TBError("Found less than -1.")
if not 10 in testArray:
    logger.TBError("Expected to find 10.")


clippedArray = clipRange(testArray,maxVal=1)
for i in range(len(clippedArray)):
    if clippedArray[i] > 1:
      logger.TBError("Found greater than 1.")
if not -10 in testArray:
    logger.TBError("Expected to find -10.")


logger.TBPass("All tests passed.")
