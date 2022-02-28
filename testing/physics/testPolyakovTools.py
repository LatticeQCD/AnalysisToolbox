# 
# testPolyakovTools.py
# 
# D. Clarke
# 
# Make sure Polyakov loop observables are calculated correctly. This test assumes that what was in the AnalysisToolbox
# as of 28 Feb 2021 was correct.
#

from latqcdtools.physics.polyakovTools import *
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.base.check import print_results

EPSILON=1e-15
NBLOCKS=40

Ns = 32
Nt = 8

pt = polyakovTools(Ns, Nt)

stream, conf, ReParr, ImParr = np.loadtxt('polyakovTable.d',unpack=True)

absPm    , absPe     = jackknife( pt.absPLoop        , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
absP2m   , absP2e    = jackknife( pt.absPLoop2       , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
chiPm    , chiPe     = jackknife( pt.Suscept         , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
T3chiPm  , T3chiPe   = jackknife( pt.T3Suscept       , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
chiRePm  , chiRePe   = jackknife( pt.ReSuscept       , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
chiImPm  , chiImPe   = jackknife( pt.ImSuscept       , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
T3chiRePm, T3chiRePe = jackknife( pt.T3ReSuscept     , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
T3chiImPm, T3chiImPe = jackknife( pt.T3ImSuscept     , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
RePxImPm , RePxImPe  = jackknife( pt.ReTimesIm       , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
Rm       , Re        = jackknife( pt.RatSuscFunction , [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
RAm      , RAe       = jackknife( pt.RatSuscFunctionA, [ReParr, ImParr], numb_blocks=NBLOCKS, conf_axis=1 )
Re2m     , Re2e      = jackknife( pt.VP2             ,  ReParr         , numb_blocks=NBLOCKS, conf_axis=1 )

print_results(absPm    , 0.0035226075714959704 ,absPe    , 7.230320110513578e-06 , "<|P|>"      , EPSILON)
print_results(absP2m   , 1.24087118252317e-05  ,absP2e   , 5.093993657785888e-08 , "<|P|^2>"    , EPSILON)
print_results(chiPm    , 0.05859092118989421   ,chiPe    , 0.00037433242806799084, "chi_P"      , EPSILON)
print_results(T3chiPm  , 0.00011443539294901213,T3chiPe  , 7.311180235702946e-07 , "T^3 chi_P"  , EPSILON)
print_results(chiRePm  , 0.06794029308913838   ,chiRePe  , 0.00048557435515950625, "chi_ReP"    , EPSILON)
print_results(chiImPm  , 0.06128526253171901   ,chiImPe  , 0.0004937885020669289 , "chi_ImP"    , EPSILON)
print_results(T3chiRePm, 0.0001326958849397234 ,T3chiRePe, 9.483874124209106e-07 , "T^3 chi_ReP", EPSILON)
print_results(T3chiImPm, 0.00011969777838226369,T3chiImPe, 9.644306680994706e-07 , "T^3 chi_ImP", EPSILON)
print_results(RePxImPm , -0.0006513185051568043,RePxImPe , 0.0007752663347590344 , "<ReP ImP>"  , EPSILON)
print_results(Rm       , 0.902013095819781     ,Re       , 0.008395693339531076  , "R"          , EPSILON)
print_results(RAm      , 0.8623775437032031    ,RAe      , 0.003329838469183075  , "R_A"        , EPSILON)
print_results(Re2m     , 0.4039126275203667    ,Re2e     , 0.0016378670645917858 , "<ReP^2>"    , EPSILON)
