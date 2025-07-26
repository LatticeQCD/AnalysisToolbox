latqcdtools.physics.constants
=============

```Python
GeV_to_fminv(x) -> float:
'''
'''
```
```Python
GeVinv_to_fm(x) -> float:
'''
'''
```
```Python
M_K0_phys(year=2022, units='MeV', returnErr=False):
'''
Physical value of K0 mass. 
'''
```
```Python
M_Kpm_phys(year=2022, units='MeV', returnErr=False):
'''
Physical value of K0 mass. 
'''
```
```Python
M_mu_phys(year=2022, units='MeV', returnErr=False):
'''
Physical value of the muon mass. 
'''
```
```Python
M_neutron_phys(year=2024, units='MeV', returnErr=False):
'''
Physical value of the neutron mass. 
'''
```
```Python
M_pi0_phys(year=2022, units='MeV', returnErr=False):
'''
Physical value of the pi0 mass. 
'''
```
```Python
M_pipm_phys(year=2022, units='MeV', returnErr=False):
'''
Physical value of the pi+/- mass. 
'''
```
```Python
M_proton_phys(year=2024, units='MeV', returnErr=False):
'''
Physical value of the proton mass. 
'''
```
```Python
M_rho_phys(year=2022, units='MeV', returnErr=False):
'''
Physical value of the rho mass. 
'''
```
```Python
MeV_to_fminv(x) -> float:
'''
'''
```
```Python
MeVinv_to_fm(x) -> float:
'''
'''
```
```Python
Rproton_phys(year=2018, units='fm', returnErr=False):
'''
Physical value of proton charge radius. 
'''
```
```Python
_separatePrefix(units):
'''
'''
```
```Python
alpha_e(year=2018, returnErr=False):
'''
Fine structure constant. 
'''
```
```Python
convert(x, unit1, unit2) -> float:
'''
General method for doing unit conversions. He knows about scientific prefixes like G, M, and so on.
If the unit ends in 'inv', it is interpreted as 1/unit.

Args:
    x (float): measurement in [unit1]. 
    unit1 (str): Original units.
    unit2 (str): Target units.

Returns:
    float: measurement in [unit2]. 
'''
```
```Python
fk_phys(year=2019, units='MeV', returnErr=False):
'''
Physical value of Kaon decay constant, f_K+/-. Scaled by sqrt(2.), which is what HotQCD usually does. 
'''
```
```Python
fm_to_GeVinv(x) -> float:
'''
'''
```
```Python
fm_to_MeVinv(x) -> float:
'''
'''
```
```Python
fminv_to_MeV(x) -> float:
'''
'''
```
```Python
fpi_phys(year=2018, units='MeV', returnErr=False):
'''
Physical value of the pion decay constant, f_pi+/-. 
'''
```
```Python
frho_phys(year=2017, units='GeV', returnErr=False):
'''
Physical value of the rho decay constant. 
'''
```
```Python
lambda_MSbar_phys(year=2021, units='MeV', returnErr=False):
'''
Physical value of MS-bar lambda parameter. 
'''
```
```Python
r0_phys(year=2014, units='fm', returnErr=False):
'''
Physical value of Sommer scale r0. 
'''
```
```Python
r1_phys(year=2010, units='fm', returnErr=False):
'''
Physical value of Sommer scale r1. 
'''
```
```Python
w0_phys(year=2013, units='fm', returnErr=False):
'''
Gradient flow scale w0.
'''
```
```Python
class physicalConstant(name, scale, units):
'''
'''
```
