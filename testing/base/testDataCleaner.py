# 
# test_dataCleaner.py
# 
# D. Clarke
# 
# Some methods to test the cleanData.py module. 
# 


import numpy as np
from latqcdtools.base.cleanData import clipRange, excludeAtCol, restrictAtCol, intersectAtCol
from latqcdtools.math.math import print_results
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testDataCleaner():

    testArray = np.linspace(-10, 10, 1000)


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

    logger.TBPass("clipRange tests passed")


    testArray =np.array([[1,2,2,1,1,2,1,2,2,1],
                         [1,2,3,4,5,6,7,8,9,0]])
    correctResult=np.array([[1,1,1,1,1],[1,4,5,7,0]])
    print_results(restrictAtCol(testArray,0,1)[0],correctResult[0],text="restrictAtCol")
    print_results(excludeAtCol(testArray,0,2)[0],correctResult[0],text="excludeAtCol")


    table1 = np.array([[0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2],
                       [0,1,2,3,4,6,7,8,9,0,1,2,3,4,6,7,8,9,0,1,3,4,5,6,7,8,9]])
    table2 = np.array([[0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2],
                       [1,2,3,4,5,6,7,8,9,0,2,3,4,5,6,7,8,9,0,1,3,4,5,6,7,8,9]])
    for stream in [0,1,2]:
        temp1 = restrictAtCol(table1,0,stream)
        temp2 = restrictAtCol(table2,0,stream)
        res1, res2 = intersectAtCol(temp1,temp2,1) 
        print_results(res1,res2,text='intersectAtCol stream '+str(stream))


if __name__ == '__main__':
    testDataCleaner()