# 
# testReadWriteConf.py                                                               
# 
# D. Clarke
# 
# To test correct read/write of configurations.
# 

from latqcdtools.interfaces.gauge import fieldNERSC
import latqcdtools.base.logger as logger

gauge = fieldNERSC(Ns=8, Nt=4)

gauge.readConf('nersc.l8t4b3360')

logger.TBPass('All tests passed!')