# 
# testReadWriteConf.py                                                               
# 
# D. Clarke
# 
# To test correct read/write of configurations.
# 

import numpy as np
from latqcdtools.interfaces.confReader import NERSCReader, ILDGReader, DynisoReader
from latqcdtools.base.utilities import timer
from latqcdtools.math.math import rel_check
from latqcdtools.testing import concludeTest
import latqcdtools.base.logger as logger


def testReadWriteConf():

    timing = timer()

    reader = NERSCReader(Ns=8, Nt=4)
    gauge1 = reader.readConf('../../datasets/nersc.l8t4b3360')
    timing.printTiming()
    gauge1.checkSU3()
    timing.printTiming()

    reader = ILDGReader(Ns=8, Nt=4)
    gauge2 = reader.readConf('../../datasets/ildg.l8t4b3360')
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

    # The dyniso files were written by dyniso's save_lattice code from the NERSC conf, in single precision. The
    # ascii format keeps 6 significant digits.
    for fmt, prec in [('comp',1e-6), ('milc',1e-6), ('ascii',1e-5)]:
        reader = DynisoReader(Ns=8, Nt=4)
        gauge3 = reader.readConf('../../datasets/dyniso_'+fmt+'.l8t4b3360')
        timing.printTiming()
        gauge3.checkSU3()
        diff = np.max(np.abs(gauge3.field-gauge1.field))
        if diff > prec:
            logger.TBFail('dyniso',fmt,'differs from NERSC conf by',diff)
            lpass = False
        else:
            logger.TBPass('dyniso',fmt,'agrees with NERSC conf.')

    concludeTest(lpass)

if __name__ == '__main__':
    testReadWriteConf()