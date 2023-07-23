# 
# testLogger.py
# 
# D. Clarke
# 
# Test some of the methods in the logger module.
# 

import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import parallel_function_eval

logger.set_log_level('ALL')
logger.createLogFile('testLogger.log')

def testLogger(i):

    logger.debug('Never again to revisit my boyhood in Surrey!')
    logger.details('Romping with my school chums')
    logger.progress('among the fens and spinneys...')
    logger.info('til the twilight bathed the hedgerows')
    logger.warn('like a lambent flame.')

if __name__ == '__main__':
    parallel_function_eval(testLogger,range(2),nproc=2)