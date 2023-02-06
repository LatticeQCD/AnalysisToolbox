# 
# testReadWriteConf.py                                                               
# 
# D. Clarke
# 
# To test correct read/write of configurations.
# 

from latqcdtools.interfaces.confReader import NERSCReader
import latqcdtools.base.logger as logger

reader = NERSCReader(Ns=8, Nt=4)

reader.readConf('nersc.l8t4b3360')

logger.TBPass('All tests passed!')