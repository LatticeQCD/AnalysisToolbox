# 
# testPolyakovTools.py
# 
# D. Clarke
# 
# Make sure Polyakov loop observables are calculated correctly. This test assumes that what was in the AnalysisToolbox
# as of 28 Feb 2021 was correct.
#

from latqcdtools.physics.polyakovTools import polyakovTools
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.base.readWrite import readTable
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testPolyakovTools():

    lpass = True

    EPSILON=1e-15
    NBLOCKS=40

    Ns = 32
    Nt = 8

    pt = polyakovTools(Ns, Nt)

    _, _, ReParr, ImParr = readTable('polyakovTable.txt')

    absPm    , absPe     = jackknife( pt.absPLoop        , [ReParr, ImParr], numb_blocks=NBLOCKS )
    absP2m   , absP2e    = jackknife( pt.absPLoop2       , [ReParr, ImParr], numb_blocks=NBLOCKS )
    chiPm    , chiPe     = jackknife( pt.Suscept         , [ReParr, ImParr], numb_blocks=NBLOCKS )
    T3chiPm  , T3chiPe   = jackknife( pt.T3Suscept       , [ReParr, ImParr], numb_blocks=NBLOCKS )
    chiRePm  , chiRePe   = jackknife( pt.ReSuscept       , [ReParr, ImParr], numb_blocks=NBLOCKS )
    chiImPm  , chiImPe   = jackknife( pt.ImSuscept       , [ReParr, ImParr], numb_blocks=NBLOCKS )
    T3chiRePm, T3chiRePe = jackknife( pt.T3ReSuscept     , [ReParr, ImParr], numb_blocks=NBLOCKS )
    T3chiImPm, T3chiImPe = jackknife( pt.T3ImSuscept     , [ReParr, ImParr], numb_blocks=NBLOCKS )
    RePxImPm , RePxImPe  = jackknife( pt.ReTimesIm       , [ReParr, ImParr], numb_blocks=NBLOCKS )
    Rm       , Re        = jackknife( pt.RatSuscFunction , [ReParr, ImParr], numb_blocks=NBLOCKS )
    RAm      , RAe       = jackknife( pt.RatSuscFunctionA, [ReParr, ImParr], numb_blocks=NBLOCKS )
    Re2m     , Re2e      = jackknife( pt.ReP2            , [ReParr, ImParr], numb_blocks=NBLOCKS )

    lpass *= print_results( absPm    , 0.0035226075714959704 , absPe    , 7.230320110513578e-06 , "<|P|>"      , EPSILON )
    lpass *= print_results( absP2m   , 1.24087118252317e-05  , absP2e   , 5.093993657785888e-08 , "<|P|^2>"    , EPSILON )
    lpass *= print_results( chiPm    , 0.05859092118989421   , chiPe    , 0.00037433242806799084, "chi_P"      , EPSILON )
    lpass *= print_results( T3chiPm  , 0.00011443539294901213, T3chiPe  , 7.311180235702946e-07 , "T^3 chi_P"  , EPSILON )
    lpass *= print_results( chiRePm  , 0.06794029308913838   , chiRePe  , 0.00048557435515950625, "chi_ReP"    , EPSILON )
    lpass *= print_results( chiImPm  , 0.06128526253171901   , chiImPe  , 0.0004937885020669289 , "chi_ImP"    , EPSILON )
    lpass *= print_results( T3chiRePm, 0.0001326958849397234 , T3chiRePe, 9.483874124209106e-07 , "T^3 chi_ReP", EPSILON )
    lpass *= print_results( T3chiImPm, 0.00011969777838226369, T3chiImPe, 9.644306680994706e-07 , "T^3 chi_ImP", EPSILON )
    lpass *= print_results( RePxImPm , -0.0006513185051568043, RePxImPe , 0.0007752663347590344 , "<ReP ImP>"  , EPSILON )
    lpass *= print_results( Rm       , 0.902013095819781     , Re       , 0.008395693339531076  , "R"          , EPSILON )
    lpass *= print_results( RAm      , 0.8623775437032031    , RAe      , 0.003329838469183075  , "R_A"        , EPSILON )
    lpass *= print_results( Re2m     , 0.4039126275203667    , Re2e     , 0.0016378670645917858 , "<ReP^2>"    , EPSILON )

    concludeTest(lpass)


if __name__ == '__main__':
    testPolyakovTools()