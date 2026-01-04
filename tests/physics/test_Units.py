# 
# testUnits.py                                                               
# 
# D. Clarke 
# 
# Check some of the unit conversion methods. 
# 

import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import convert, fm_to_MeVinv, MeVinv_to_fm, fm_to_GeVinv, GeVinv_to_fm, \
    M_mu_phys, M_pi0_phys, M_pipm_phys, fk_phys, frho_phys, lambda_MSbar_phys, cms, sqrtG, M_e_phys, M_proton_phys
from latqcdtools.testing import print_results, concludeTest

logger.set_log_level('INFO')

def testUnits():

    lpass = True

    lpass *= print_results(convert(1,'km','m'),1000,text='[km] to [m]')
    lpass *= print_results(convert(60,'s','min'),1,text='[s] to [min]')
    lpass *= print_results(convert(3600,'s','h'),1,text='[s] to [h]')
    lpass *= print_results(convert(31556926.080000006,'s','y'),1,text='[s] to [y]')
    lpass *= print_results(convert(8765.8128,'h','y'),1,text='[h] to [y]')
    lpass *= print_results(convert(8765.8128,'h','y'),1,text='[h] to [y]')
    lpass *= print_results(convert(94e6,'mi','m'),1.51278336e11,text='[mi] to [m]')
    lpass *= print_results(convert(cms,'m/s','mi/h'),670616629.3843951,text='[m/s] to [mi/h]')
    lpass *= print_results(convert(10632,'kWh/y','W'),1212.8938003330395,text='[kWh/y] to [W]')
    lpass *= print_results(convert(1,'Wh','BTU'),3.412,text='[Wh] to [BTU]')
    lpass *= print_results(convert(1,'Wh','J'),3600,text='[Wh] to [J]')
    lpass *= print_results(convert(1,'eV','J'),1.602176634e-19,text='[eV] to [J]')
    lpass *= print_results(convert(1,'eV','K'),11604.518121550082,text='[eV] to [K]')
    lpass *= print_results(convert(10,'kK','degC'),9726.85,text='[K] to [degC]')
    lpass *= print_results(convert(12,'kdegC','degF'),21632,text='[degC] to [degF]')

    # Identity operation checks
    lpass *= print_results( convert( convert(1,'s'   ,'min' ),'min' ,'s'   ), 1, text='identity: [s], [min]' )
    lpass *= print_results( convert( convert(1,'s'   ,'h'   ),'h'   ,'s'   ), 1, text='identity: [s], [h]' )
    lpass *= print_results( convert( convert(1,'s'   ,'y'   ),'y'   ,'s'   ), 1, text='identity: [s], [y]' )
    lpass *= print_results( convert( convert(1,'h'   ,'y'   ),'y'   ,'h'   ), 1, text='identity: [h], [y]' )
    lpass *= print_results( convert( convert(1,'mi'  ,'m'   ),'m'   ,'mi'  ), 1, text='identity: [m], [mi]' )
    lpass *= print_results( convert( convert(1,'mi'  ,'ft'  ),'ft'  ,'mi'  ), 1, text='identity: [mi], [ft]' )
    lpass *= print_results( convert( convert(1,'ft'  ,'m'   ),'m'   ,'ft'  ), 1, text='identity: [ft], [m]' )
    lpass *= print_results( convert( convert(1,'m/s' ,'mi/h'),'mi/h','m/s' ), 1, text='identity: [mi/h], [m/s]' )
    lpass *= print_results( convert( convert(1,'Wh/y','W'   ),'W'   ,'Wh/y'), 1, text='identity: [Wh/y], [W]' )
    lpass *= print_results( convert( convert(1,'Wh'  ,'BTU' ),'BTU' ,'Wh'  ), 1, text='identity: [Wh], [BTU]' )
    lpass *= print_results( convert( convert(1,'Wh'  ,'J'   ),'J'   ,'Wh'  ), 1, text='identity: [Wh], [J]' )
    lpass *= print_results( convert( convert(1,'eV'  ,'J'   ),'J'   ,'eV'  ), 1, text='identity: [eV], [J]' )
    lpass *= print_results( convert( convert(1,'eV'  ,'K'   ),'K'   ,'eV'  ), 1, text='identity: [eV], [K]' )
    lpass *= print_results( convert( convert(1,'eV'  ,'minv'),'minv','eV'  ), 1, text='identity: [eV], [1/m]' )
    lpass *= print_results( convert( convert(1,'K'   ,'degC'),'degC','K'   ), 1, text='identity: [K], [degC]' )
    lpass *= print_results( convert( convert(1,'K'   ,'degF'),'degF','K'   ), 1, text='identity: [K], [degF]' )
    lpass *= print_results( convert( convert(1,'degF','degC'),'degC','degF'), 1, text='identity: [degF], [degC]' )
    lpass *= print_results( fm_to_MeVinv( MeVinv_to_fm(1) ), 1, text='identity: [fm], [1/MeV]' )
    lpass *= print_results( fm_to_GeVinv( GeVinv_to_fm(1) ), 1, text='identity: [fm], [1/GeV]' )

    # Tests of some constants
    lpass *= print_results(M_mu_phys(2020,"MeV"), 105.6583745,text="m_mu [MeV]")
    lpass *= print_results(M_mu_phys(2020,"GeV"), 0.1056583745,text="m_mu [GeV]")
    lpass *= print_results(M_mu_phys(2020,"fminv"), 0.5354481943753351,text="m_mu [1/fm]")
    lpass *= print_results(M_pi0_phys(2022,"MeV"), 134.9768,text="m_pi0 [MeV]")
    lpass *= print_results(M_pi0_phys(2022,"GeV"), 0.1349768,text="m_pi0 [GeV]")
    lpass *= print_results(M_pi0_phys(2022,"fminv"), 0.6840260810804044, text="m_pi0 [1/fm]")
    lpass *= print_results(M_pipm_phys(2022,"MeV"), 139.57039,text="m_pipm [MeV]")
    lpass *= print_results(M_pipm_phys(2022,"GeV"), 0.13957039000000002, text="m_pipm [GeV]")
    lpass *= print_results(fk_phys(2019,"GeV"), 0.11009652583074543,text="fk [GeV]")
    lpass *= print_results(fk_phys(2019,"MeV"), 110.09652583074543,text="fk [MeV]")
    lpass *= print_results(fk_phys(2019,"fminv"), 0.5579395503862317,text="fk [1/fm]")
    lpass *= print_results(frho_phys(2017,"GeV"), 0.21,text="frho [GeV]")
    lpass *= print_results(frho_phys(2017,"MeV"), 210.0,text="frho [MeV]")
    lpass *= print_results(frho_phys(2017,"fminv"), 1.06422345934179,text="frho [1/fm]")
    lpass *= print_results(lambda_MSbar_phys(2021,"GeV"), 0.339,text="lambda_MSbar [GeV]")
    lpass *= print_results(lambda_MSbar_phys(2021,"MeV"), 339,text="lambda_MSbar [MeV]")
    lpass *= print_results(lambda_MSbar_phys(2021,"fminv"), 1.7179607272231752,text="lambda_MSbar [1/fm]")

    # Compute G_N two ways
    Gm = 6.67430e-11    # Units: m^3 kg^-1 s^-2
    for i in range(4):
        Gm = convert(Gm,"m","s") 
    Gm = convert( convert(Gm,"m","GeVinv"), "Jinv", "GeVinv" )
    lpass *= print_results(sqrtG(2024,"GeVinv")**2, Gm,text="G_B [1/GeV^2]")

    # Some more tests of natural unit conversions. Test against reference https://public.websites.umich.edu/~jwells/Scholardox/A3.pdf.
    # They rounded to 3 significant figures, so that's the best accuracy I can get. 
    lpass *= print_results(convert(M_e_phys(units="MeV")     ,"MeV","fminv"),1/386  ,text="m_e natural units",prec=1e-3)
    lpass *= print_results(convert(M_proton_phys(units="MeV"),"MeV","fminv"),1/0.210,text="m_p natural units",prec=3e-3)

    # More tests of natural unit conversions. Test against Bastian Brandt's lecture notes.
    lpass *= print_results(convert(1,'h','eVinv'), 5.46936285835437e+18,text='[h] to [1/eV]')
    lpass *= print_results(convert(1,'eVinv','s'), 6.582119514160693e-16,text='[1/eV] to [s]')
    lpass *= print_results(convert(1,'kg','eV'), 5.609588603804453e+35, text='[kg] to [eV]')

    # Some natural unit identity operation checks
    lpass *= print_results(convert(convert(1,'eVinv','s'),'s','eVinv'),1,text='identity: [1/eV], [s]')
    lpass *= print_results(convert(convert(1,'eVinv','h'),'h','eVinv'),1,text='identity: [1/eV], [h]')
    lpass *= print_results(convert(convert(1,'h'    ,'m'),'m','h'    ),1,text='identity: [h], [m]'   )
    lpass *= print_results(convert(convert(1,'eV'   ,'g'),'g','eV'   ),1,text='identity: [eV], [g]'  )

    # Also a sanity check for the prefixes and physical parameter class
    for prefix in ["Q","R","Y","Z","E","P","T","G","M","k","h"]:
        lpass *= print_results(convert(M_proton_phys(units=f"{prefix}eV"),f"{prefix}eV","fminv"),1/0.210,
                               text=f"m_p natural, prefix={prefix}",prec=3e-3)

    concludeTest(lpass)


if __name__ == '__main__':
    testUnits()