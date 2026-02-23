# 
# testScales.py                                                               
# 
# D. Clarke
# 
# Test the lattice QCD reference scales. 
#

import numpy as np
from latqcdtools.physics.referenceScales import a_times_fk, a_div_r1, r0_div_a, sqrtt0_div_a
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.physics.constants import M_mu_phys, M_K0_phys, M_Kpm_phys, M_pi0_phys, M_pipm_phys,\
    M_rho_phys, frho_phys, fk_phys, fpi_phys, w0_phys, r0_phys, r1_phys, alpha_e, lambda_MSbar_phys,\
    Rproton_phys, Tpc_chiral
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')
EPSILON=1e-12


def testScales():

    lpass = True

    # Test different parameterization years. This assumes that what was in the AnalysisToolbox 
    # as of 28 Feb 2021 was correct.
    lpass *= print_results( a_times_fk(6.500, 2021), 0.07826294754259573, text="a fK 2021" , prec=EPSILON )
    lpass *= print_results( a_times_fk(6.500, 2014), 0.07870391862551822, text="a fK 2014" , prec=EPSILON )
    lpass *= print_results( a_times_fk(6.500, 2012), 0.07855403707936452, text="a fK 2012" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2021)  , 0.17296284805472273, text="a/r1 2021" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2018)  , 0.17293304059842668, text="a/r1 2018" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2014)  , 0.17282128740434255, text="a/r1 2014" , prec=EPSILON )
    lpass *= print_results( a_div_r1(7.500, 2012)  , 0.17512946173066155, text="a/r1 2012" , prec=EPSILON )
    lpass *= print_results( r0_div_a(6.800, 2017)  , 16.440532234860257 , text="r0/a"      , prec=EPSILON )
    lpass *= print_results( sqrtt0_div_a(6.800)    , 5.524714874614185  , text="sqrt(t0)/a", prec=EPSILON )

    # Test some physical scales. This should be improved to be more careful in the future, but as is,
    # it's better than nothing.
    lpass *= print_results( 105.6583755      , M_mu_phys(year=2022,units="MeV")        , text="muon mass 2022"    )
    lpass *= print_results( 134.9768         , M_pi0_phys(year=2022,units="MeV")       , text="pi0 mass 2022"     )
    lpass *= print_results( 139.57039        , M_pipm_phys(year=2022,units="MeV")      , text="pi+/- mass 2022"   )
    lpass *= print_results( 497.611          , M_K0_phys(year=2022,units="MeV")        , text="K0 mass 2022"      )
    lpass *= print_results( 493.677          , M_Kpm_phys(year=2022,units="MeV")       , text="K+/- mass 2022"    )
    lpass *= print_results( 775.26           , M_rho_phys(year=2022,units="MeV")       , text="rho mass 2022"     )
    lpass *= print_results(   0.21           , frho_phys(year=2017,units="GeV")        , text="f_rho 2017"        )
    lpass *= print_results( 155.7/np.sqrt(2.), fk_phys(year=2019,units="MeV")          , text="f_K 2019"          )
    lpass *= print_results( 130.50           , fpi_phys(year=2018,units="MeV")         , text="f_pi 2018"         )
    lpass *= print_results(   0.1715         , w0_phys(year=2013,units="fm")           , text="w0 2018"           )
    lpass *= print_results(   0.3106         , r1_phys(year=2010,units="fm")           , text="r1 2010"           )
    lpass *= print_results(   0.46875752     , r0_phys(year=2014,units="fm")           , text="r0 2014"           )
    lpass *= print_results(   7.2973525693e-3, alpha_e(year=2018)                      , text="alpha_e 2018"      )
    lpass *= print_results( 339              , lambda_MSbar_phys(year=2021,units="MeV"), text="lambda_MSbar 2021" )
    lpass *= print_results(   0.8414         , Rproton_phys(year=2018,units="fm")      , text="R_proton 2021"     )
    lpass *= print_results( 157.2555509      , Tpc_chiral(year=2020,units="MeV")       , text="Tpc_chiral 2020"   )

    concludeTest(lpass)

if __name__ == '__main__':
    testScales()
