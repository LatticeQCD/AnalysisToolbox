# Reference scales and physical constants

## Scale setting

In the module 
```Python
latqcdtools.physics.referenceScales
```
you find lots of functions for scale setting, along with references for each scale determination. 
Supported scales include: $f_K$, $t_0$, $r_0$, $r_1$, and $r_1m_s$. 

You can get scales in lattice units
through methods like `a_times_fk(beta: float, year)`, which takes the bare coupling 
as input and returns $a\,f_K(\beta)$
using a parameterization specified by the year `year`. In this particular instance, there
exists parameterizations from 2021, 2014, and 2012.
Please note that the parameterizations for $a$ have closed intervals $[\beta_{\rm min},~\beta_{\rm max}]$ in
which they are valid.


## Physical constants

In the module
```Python
latqcdtools.physics.constants
```
one finds many constants of nature. In general these have the form
```Python
M_rho_phys(year=2022,units="MeV",returnErr=False)
```
The default year is the most recent measurement the authors have implemented. This number is usually
taken from the PDG, but it may be taken from other papers. We always try to list the source of the measurement.
If `returnErr=True`, this will return a tuple of the form `mean, error`. The list of included physical
constants expands as our authors work on more and more projects.

What you might also find helpful in this module is a method for conversion between many kinds of units.
This is carried out with
```Python
convert(x,unit1,unit2)
```
which converts the number `x` from `unit1` to units `unit2`. Note that `unit1` and `unit2` know about
scientific prefixes. This includes some of the standard conversions in high energy physics, like
fm to MeV$^{-1}$. But there are also time units in case you work on something touching cosmology
and energy units in case you want to optimize your electricity bill.

