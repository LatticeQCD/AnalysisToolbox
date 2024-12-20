latqcdtools.physics.constants
=============

`GeV_to_fminv(x) -> float`


`GeVinv_to_fm(x) -> float`


`M_K0_phys(year=2022, units='MeV', returnErr=False)`
 
    Physical value of K0 mass. 
    
`M_Kpm_phys(year=2022, units='MeV', returnErr=False)`
 
    Physical value of K0 mass. 
    
`M_mu_phys(year=2022, units='MeV', returnErr=False)`
 
    Physical value of the muon mass. 
    
`M_pi0_phys(year=2022, units='MeV', returnErr=False)`
 
    Physical value of the pi0 mass. 
    
`M_pipm_phys(year=2022, units='MeV', returnErr=False)`
 
    Physical value of the pi+/- mass. 
    
`M_rho_phys(year=2022, units='MeV', returnErr=False)`
 
    Physical value of the rho mass. 
    
`MeV_to_fminv(x) -> float`


`MeVinv_to_fm(x) -> float`


`Rproton_phys(year=2018, units='fm', returnErr=False)`
 
    Physical value of proton charge radius. 
    
`_separatePrefix(units)`


`alpha_e(year=2018, returnErr=False)`
 
    Fine structure constant. 
    
`convert(x, unit1, unit2) -> float`
 
    General method for doing unit conversions. He knows about scientific prefixes like G, M, and so on.
    If the unit ends in 'inv', it is interpreted as 1/unit.

    Args:
        x (float): measurement in [unit1]. 
        unit1 (str): Original units.
        unit2 (str): Target units.

    Returns:
        float: measurement in [unit2]. 
    
`fk_phys(year=2019, units='MeV', returnErr=False)`
 
    Physical value of Kaon decay constant, f_K+/-. Scaled by sqrt(2.), which is what HotQCD usually does. 
    
`fm_to_GeVinv(x) -> float`


`fm_to_MeVinv(x) -> float`


`fpi_phys(year=2018, units='MeV', returnErr=False)`

    Physical value of the pion decay constant, f_pi+/-. 
    
`frho_phys(year=2017, units='GeV', returnErr=False)`
 
    Physical value of the rho decay constant. 
    
`lambda_MSbar_phys(year=2021, units='MeV', returnErr=False)`
 
    Physical value of MS-bar lambda parameter. 
    
`r0_phys(year=2014, units='fm', returnErr=False)`
 
    Physical value of Sommer scale r0. 
    
`r1_phys(year=2010, units='fm', returnErr=False)`
 
    Physical value of Sommer scale r1. 
    
`w0_phys(year=2013, units='fm', returnErr=False)`
 
    Gradient flow scale w0.
    
`physicalConstant(name, scale, units)`


