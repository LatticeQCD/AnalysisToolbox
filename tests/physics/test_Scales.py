# 
# testScales.py                                                               
# 
# D. Clarke
# 
# Test the lattice QCD reference scales. This test assumes that what was in the AnalysisToolbox as of 28 Feb 2021
# was correct.
#

from latqcdtools.physics.referenceScales import a_times_fk, a_div_r1, r0_div_a, sqrtt0_div_a
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.physics.constants import M_mu_phys
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')
EPSILON=1e-15


def testScales():

    lpass = True

    # In each case the physical scale is kept fixed and the parameterization year is varied.
    lpass *= print_results( a_times_fk(6.500, 2021), 0.07826294754259573, text="a fK 2021" , prec=EPSILON )
    lpass *= print_results( a_times_fk(6.500, 2014), 0.07870391862551822, text="a fK 2014" , prec=EPSILON )
    lpass *= print_results( a_times_fk(6.500, 2012), 0.07855403707936452, text="a fK 2012" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2021)  , 0.17296284805472273, text="a/r1 2021" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2018)  , 0.17293304059842668, text="a/r1 2018" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2014)  , 0.17282128740434255, text="a/r1 2014" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2012)  , 0.17512946173066155, text="a/r1 2012" , prec=EPSILON )
    lpass *= print_results( r0_div_a(6.800, 2017)  , 16.440532234860257 , text="r0/a"      , prec=EPSILON )
    lpass *= print_results( sqrtt0_div_a(6.800)    , 5.524714874614185  , text="sqrt(t0)/a", prec=EPSILON )

    lpass *= print_results( 105.6583755, M_mu_phys(year=2022,units="MeV"), text="muon mass 2022" )

    concludeTest(lpass)

if __name__ == '__main__':
    testScales()