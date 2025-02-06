# 
# testStatPhys.py                                                               
# 
# D. Clarke
# 
# Test some of the statistical physics methods.
#


from latqcdtools.physics.statisticalPhysics import O2_3d, O3_3d, O4_3d, Z2_3d, Z2_2d, S3_2d, S4_2d, reweight
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable
from latqcdtools.testing import concludeTest, gaudif_results
from latqcdtools.statistics.statistics import std_mean
from latqcdtools.statistics.jackknife import jackknife
import glob
import numpy as np


logger.set_log_level('INFO')


def testStatPhys():

    lpass = True

    univ = O2_3d()
    lpass *= univ.hyperscalingCheck()
    univ = O3_3d()
    lpass *= univ.hyperscalingCheck(tol=1e-2)
    univ = O4_3d()
    lpass *= univ.hyperscalingCheck()
    univ = Z2_3d()
    lpass *= univ.hyperscalingCheck()
    univ = Z2_2d()
    lpass *= univ.hyperscalingCheck()
    if univ.Tc() != 2/np.log(1+np.sqrt(2)):
        logger.TBFail('Z2_2d Tc')
    univ = S3_2d()
    lpass *= univ.hyperscalingCheck()
    univ = S4_2d()
    lpass *= univ.hyperscalingCheck()

    # A test of the reweighter where we check against 16x16 Ising model data

    V        = 16**2
    p0       = 1/2.3
    M0, act0 = readTable('RWcontrol/T2.3.d')
    S        = act0*V/2
    
    def RW(data,xRW,x0):
        """ Change the order of arguments to be compatible with jackknife. """
        X = data[0]
        Y = data[1]
        return reweight(X,xRW,x0,Y)
    
    for filename in glob.iglob('RWcontrol/T*.d'):
    
        T          = float(filename[11:-2])
        pRW        = 1/T
        M, _       = readTable(filename)
        MRWm, MRWe = jackknife( RW, [M0, S], args=(pRW,p0) )
        Mm  , Me   = jackknife( std_mean, M )
        lpass     *= gaudif_results(MRWm,MRWe,Mm,Me,text='T='+str(T)) 

    concludeTest(lpass)


if __name__ == '__main__':
    testStatPhys()