latqcdtools.physics.constants
=============

```python
GeV_to_fminv(x) -> float
```
```python
GeVinv_to_fm(x) -> float
```
```python
M_K0_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_Kpm_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_b_phys(year=2024, units='GeV', returnErr=False, world='nature')
```
```python
M_c_phys(year=2024, units='GeV', returnErr=False, world='nature')
```
```python
M_d_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```python
M_e_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```python
M_mu_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_neutron_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```python
M_phi_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_pi0_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_pipm_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_proton_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```python
M_rho_phys(year=2022, units='MeV', returnErr=False, world='nature')
```
```python
M_s_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```python
M_t_phys(year=2024, units='GeV', returnErr=False, world='nature')
```
```python
M_u_phys(year=2024, units='MeV', returnErr=False, world='nature')
```
```python
MeV_to_fminv(x) -> float
```
```python
MeVinv_to_fm(x) -> float
```
```python
Rproton_phys(year=2018, units='fm', returnErr=False, world='nature'):
"""
Physical value of proton charge radius. 
"""
```
```python
Tpc_chiral(year=2020, units='MeV', returnErr=False, world='Nf21'):
"""
Chiral crossover temperature at zero net-baryon chemical potential, as estimated using the
chiral susceptibility. This differs from other methods, as different years combine results
from multiple groups. 
"""
```
```python
_separatePrefix(units)
```
```python
alpha_e(year=2018, returnErr=False, world='nature'):
"""
Fine structure constant. 
"""
```
```python
convert(x, unit1, unit2) -> float:
"""
General method for doing unit conversions. He knows about scientific prefixes like G, M, and so on.
If the unit ends in 'inv', it is interpreted as 1/unit. He also knows about natural units. You can
only convert away one power or inverse power of a unit at a time.

Args:
    x (float): measurement in [unit1]. 
    unit1 (str): Original units.
    unit2 (str): Target units.

Returns:
    float: measurement in [unit2]. 
"""
```
```python
fk_phys(year=2019, units='MeV', returnErr=False, world='nature'):
"""
Physical value of Kaon decay constant, f_K+/-. Scale by sqrt(2.). 
"""
```
```python
fm_to_GeVinv(x) -> float
```
```python
fm_to_MeVinv(x) -> float
```
```python
fminv_to_MeV(x) -> float
```
```python
fphi_phys(year=2021, units='MeV', returnErr=False, world='Nf21'):
"""
Physical value of the phi decay constant.
"""
```
```python
fpi_phys(year=2018, units='MeV', returnErr=False, world='nature'):
"""
Physical value of the pion decay constant, f_pi+/-. Scale by sqrt(2.) 
"""
```
```python
frho_phys(year=2017, units='GeV', returnErr=False, world='nature'):
"""
Physical value of the rho decay constant. 
"""
```
```python
lambda_MSbar_phys(year=2021, units='MeV', returnErr=False, world='nature'):
"""
Physical value of MS-bar lambda parameter. 
"""
```
```python
r0_phys(year=2014, units='fm', returnErr=False, world='Nf21'):
"""
Physical value of Sommer scale r0. 
"""
```
```python
r1_phys(year=2010, units='fm', returnErr=False, world='Nf21'):
"""
Physical value of Sommer scale r1. 
"""
```
```python
sqrtG(year=2024, units='GeVinv', returnErr=False, world='nature'):
"""
Square root of Newton's gravitational constant.
"""
```
```python
sqrtt0_phys(year=2017, units='fm', returnErr=False, world='Nf21'):
"""
Gradient flow scale sqrt(t0).
"""
```
```python
w0_phys(year=2013, units='fm', returnErr=False, world='Nf211'):
"""
Gradient flow scale w0.
"""
```
```python
class physicalConstant(name, scale, units):
```
