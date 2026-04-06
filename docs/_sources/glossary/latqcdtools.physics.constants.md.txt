latqcdtools.physics.constants
=============

```Python
GeV_to_fminv(x) -> float
```
```Python
GeVinv_to_fm(x) -> float
```
```Python
M_K0_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_Kpm_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_b_phys(year=2024, units='GeV', returnErr=False, world='nature')
```
```Python
M_c_phys(year=2024, units='GeV', returnErr=False, world='nature')
```
```Python
M_d_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```Python
M_e_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```Python
M_mu_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_neutron_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```Python
M_phi_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_pi0_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_pipm_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_proton_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```Python
M_rho_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```Python
M_s_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```Python
M_t_phys(year=2024, units='GeV', returnErr=False, world='nature')
```
```Python
M_u_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```Python
MeV_to_fminv(x) -> float
```
```Python
MeVinv_to_fm(x) -> float
```
```Python
Rproton_phys(year=2018, units='fm', returnErr=False, world='nature'):
'''
Physical value of proton charge radius. 
'''
```
```Python
Tpc_chiral(year=2020, units='MeV', returnErr=False, world='Nf21'):
'''
Chiral crossover temperature at zero net-baryon chemical potential, as estimated using the
chiral susceptibility. This differs from other methods, as different years combine results
from multiple groups. 
'''
```
```Python
_separatePrefix(units)
```
```Python
alpha_e(year=2018, returnErr=False, world='nature'):
'''
Fine structure constant. 
'''
```
```Python
convert(x, unit1, unit2) -> float:
'''
General method for doing unit conversions. He knows about scientific prefixes like G, M, and so on.
If the unit ends in 'inv', it is interpreted as 1/unit. He also knows about natural units. You can
only convert away one power or inverse power of a unit at a time.

Args:
    x (float): measurement in [unit1]. 
    unit1 (str): Original units.
    unit2 (str): Target units.

Returns:
    float: measurement in [unit2]. 
'''
```
```Python
fk_phys(year=2019, units='MeV', returnErr=False, world='nature'):
'''
Physical value of Kaon decay constant, f_K+/-. scale by sqrt(2.). 
'''
```
```Python
fm_to_GeVinv(x) -> float
```
```Python
fm_to_MeVinv(x) -> float
```
```Python
fminv_to_MeV(x) -> float
```
```Python
fphi_phys(year=2021, units='MeV', returnErr=False, world='Nf21'):
'''
Physical value of the phi decay constant.
'''
```
```Python
fpi_phys(year=2018, units='MeV', returnErr=False, world='nature'):
'''
Physical value of the pion decay constant, f_pi+/-. 
'''
```
```Python
frho_phys(year=2017, units='GeV', returnErr=False, world='nature'):
'''
Physical value of the rho decay constant. 
'''
```
```Python
lambda_MSbar_phys(year=2021, units='MeV', returnErr=False, world='nature'):
'''
Physical value of MS-bar lambda parameter. 
'''
```
```Python
r0_phys(year=2014, units='fm', returnErr=False, world='Nf21'):
'''
Physical value of Sommer scale r0. 
'''
```
```Python
r1_phys(year=2010, units='fm', returnErr=False, world='Nf21'):
'''
Physical value of Sommer scale r1. 
'''
```
```Python
sqrtG(year=2024, units='GeVinv', returnErr=False, world='nature'):
'''
Square root of Newton's gravitational constant.
'''
```
```Python
sqrtt0_phys(year=2017, units='fm', returnErr=False, world='Nf21'):
'''
Gradient flow scale sqrt(t0).
'''
```
```Python
w0_phys(year=2013, units='fm', returnErr=False, world='Nf211'):
'''
Gradient flow scale w0.
'''
```
```Python
class physicalConstant(name, scale, units):
```
