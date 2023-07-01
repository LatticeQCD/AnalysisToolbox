# 
# testUnits.py                                                               
# 
# D. Clarke 
# 
# Check some of the unit conversion methods. 
# 

import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import convert, fm_to_MeVinv, MeVinv_to_fm, fm_to_GeVinv, GeVinv_to_fm, \
    M_mu_phys, M_pi0_phys, M_pipm_phys, fk_phys, frho_phys, lambda_MSbar_phys, separatePrefix
from latqcdtools.math.math import print_results

logger.set_log_level('INFO')

def testUnits():

    print_results(convert(1,'km','m'),1000,text='[km] to [m]')
    print_results(convert(60,'s','min'),1,text='[s] to [min]')
    print_results(convert(3600,'s','h'),1,text='[s] to [h]')
    print_results(convert(31556926.080000006,'s','y'),1,text='[s] to [y]')
    print_results(convert(8765.8128,'h','y'),1,text='[h] to [y]')
    print_results(convert(8765.8128,'W','Wh/y'),1,text='[W] to [Wh/y]')

    # Identity operation checks
    print_results( convert( convert(1,'s','min'),'min','s' ), 1, text='identity: [s], [min]' )
    print_results( convert( convert(1,'s','h'  ),'h'  ,'s' ), 1, text='identity: [s], [h]' )
    print_results( convert( convert(1,'s','y'  ),'y'  ,'s' ), 1, text='identity: [s], [y]' )
    print_results( convert( convert(1,'h','y'  ),'y'  ,'h' ), 1, text='identity: [h], [y]' )
    print_results( fm_to_MeVinv( MeVinv_to_fm(1) ), 1, text='identity: [fm], [1/MeV]' )
    print_results( fm_to_GeVinv( GeVinv_to_fm(1) ), 1, text='identity: [fm], [1/GeV]' )

    # Tests of some constants
    print_results(M_mu_phys(2020,"MeV"), 105.6583745,text="m_mu [MeV]")
    print_results(M_mu_phys(2020,"GeV"), 0.1056583745,text="m_mu [GeV]")
    print_results(M_mu_phys(2020,"fminv"), 0.5354481943753351,text="m_mu [1/fm]")
    print_results(M_pi0_phys(2022,"MeV"), 134.9768,text="m_pi0 [MeV]")
    print_results(M_pi0_phys(2022,"GeV"), 0.1349768,text="m_pi0 [GeV]")
    print_results(M_pi0_phys(2022,"fminv"), 0.6840260810804044, text="m_pi0 [1/fm]")
    print_results(M_pipm_phys(2022,"MeV"), 139.57039,text="m_pipm [MeV]")
    print_results(M_pipm_phys(2022,"GeV"), 0.13957039000000002, text="m_pipm [GeV]")
    print_results(fk_phys(2019,"GeV"), 0.11009652583074543,text="fk [GeV]")
    print_results(fk_phys(2019,"MeV"), 110.09652583074543,text="fk [MeV]")
    print_results(fk_phys(2019,"fminv"), 0.5579395503862317,text="fk [1/fm]")
    print_results(frho_phys(2017,"GeV"), 0.21,text="frho [GeV]")
    print_results(frho_phys(2017,"MeV"), 210.0,text="frho [MeV]")
    print_results(frho_phys(2017,"fminv"), 1.06422345934179,text="frho [1/fm]")
    print_results(lambda_MSbar_phys(2021,"GeV"), 0.339,text="lambda_MSbar [GeV]")
    print_results(lambda_MSbar_phys(2021,"MeV"), 339,text="lambda_MSbar [MeV]")
    print_results(lambda_MSbar_phys(2021,"fminv"), 1.7179607272231752,text="lambda_MSbar [1/fm]")


if __name__ == '__main__':
    testUnits()