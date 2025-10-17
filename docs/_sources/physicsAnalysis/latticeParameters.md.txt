# Handling Lattice Parameters

In lattice calculations we are always analyzing configurations, which depends on $m_l$, $m_s$, 
$\beta$, $N_\sigma$, and $N_\tau$. Usually the parameter combination is succinctly incorporated 
into the parameter name as, e.g. `l328f21b6390m00257m0694`. 

## Lattice parameters

The `latticeParams` class 
```Python
import latqcdtools.physics.lattice_params
``` 
collects all this information into one object, and contains some useful 
methods, for example the ability to convert one of these parameter strings into its corresponding 
float value. To instantiate this object, use
```Python
lp = latticeParams(Ns, Nt, beta, mass_l='', mass_s='',scaleType='fk', paramYear=2021, Nf='21')
```
Here `beta` and the quark masses are input strings coming from `l328f21b6390m00257m0694`. You can 
then, for instance get the float value corresponding to the string `6390` using
```Python
lp.getBeta()
```
The `scaleType` argument lets you choose what reference scale to use. How these scales vary in lattice units
as a function of $\beta$ has been calculated by the HotQCD collaboration; the `paramYear` option lets you
pick the year this function was computed. You can specify the number of flavors with `Nf`.


There is also a method to get the `l328f21b6390m00257m0694` string, `getcparams`. You can also get 
the lattice spacing in [fm] and temperature in [MeV] using `geta` and `getT`, respectively. 
Finally `paramSummary()` prints a nice summary of all parameters going into your calculation 
to the screen.

