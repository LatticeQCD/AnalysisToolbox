# 
# testReadWriteConf.py                                                               
# 
# D. Clarke
# 
# To test correct read/write of configurations.
# 

from latqcdtools.interfaces.confReader import NERSCReader
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

def testReadWriteConf():

    reader = NERSCReader(Ns=8, Nt=4)

    gauge = reader.readConf('nersc.l8t4b3360')

    print(gauge.getLink(0,0,1,1,0)) # Get the link at site (0,0,1,1) pointing in the 0 direction.

    logger.TBPass('All tests passed!')

if __name__ == '__main__':
    testReadWriteConf()