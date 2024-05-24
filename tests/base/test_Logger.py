# 
# testLogger.py
# 
# D. Clarke
# 
# Test some of the methods in the logger module.
# 

import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import parallel_function_eval

logger.createLogFile('testLogger.log')

class testClass:
    def __init__(self):
        logger.set_log_level('ALL')
        logger.debug('Never again to revisit my boyhood in Surrey!')
        logger.details('Romping with my school chums')
        logger.progress('among the fens and spinneys...')
        logger.info('til the twilight bathed the hedgerows')
        logger.warn('like a lambent flame.')
        logger.set_log_level('INFO')
    def __repr__(self) -> str:
        return "testClass"

def makeLogger(i):
    test = testClass()

def testLogger():
    parallel_function_eval(makeLogger,range(2),nproc=2)

if __name__ == '__main__':
    testLogger()
