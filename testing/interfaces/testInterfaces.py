# 
# testInterfaces.py                                                               
# 
# D. Clarke
# 
# Testing of some generic interfacing tools. 
# 

from latqcdtools.interfaces.interfaces import latexTable, redmineTable 
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

def testInterfaces():
    test = latexTable() 
    test.append([1,2,3])
    test.append([1,2,3])
    test.append([1,2,3])
    test.outputTable('testInterface.tex')
    test = redmineTable() 
    test.append([1,2,3])
    test.append([1,2,3])
    test.append([1,2,3])
    test.outputTable('testInterface.redmine')
    logger.TBPass('No problems encountered.')

if __name__ == '__main__':
    testInterfaces()