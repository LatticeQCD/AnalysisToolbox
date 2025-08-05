# 
# testReadWriteConf.py                                                               
# 
# D. Clarke
# 
# To test correct read/write of configurations.
# 

from latqcdtools.interfaces.confReader import NERSCReader, ILDGReader
from latqcdtools.base.utilities import timer
from latqcdtools.math.math import rel_check
from latqcdtools.testing import concludeTest
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

def testReadWriteConf():

    timing = timer()

    reader = NERSCReader(Ns=8, Nt=4)
    gauge1 = reader.readConf('nersc.l8t4b3360')
    timing.printTiming()
    gauge1.checkSU3()
    timing.printTiming()

    reader = ILDGReader(Ns=8, Nt=4)
    gauge2 = reader.readConf('ildg.l8t4b3360')
    timing.printTiming()
    gauge2.checkSU3()
    timing.printTiming()

    lpass = True
    for mu in range(4):
        for t in range(4):
            for z in range(8):
                for y in range(8):
                    for x in range(8):
                        link1 = gauge1.getLink(x,y,z,t,mu)
                        link2 = gauge2.getLink(x,y,z,t,mu)
                        lpass *= rel_check(link1, link2)

    concludeTest(lpass) 

if __name__ == '__main__':
    testReadWriteConf()