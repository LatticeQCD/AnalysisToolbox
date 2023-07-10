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

def testLogger(logFile='Toolbox.log'):

    logger.createLogFile(filename=logFile)
    logger.debug('Never again to revisit my boyhood in Surrey!')
    logger.details('Romping with my school chums')
    logger.progress('among the fens and spinneys...')
    logger.info('til the twilight bathed the hedgerows')
    logger.warn('like a lambent flame.')
    logger.closeLogFile()
    logger.info('This should not appear in log file.')

if __name__ == '__main__':
    testLogger()
    myFiles = ['Toolbox_parallel.log','Toolbox_parallel.log']
    parallel_function_eval(testLogger,myFiles,2)