# Handling Lattice Parameters

In lattice calculations we are always analyzing configurations, which depends on $m_l$, $m_s$, 
$\beta$, $N_\sigma$, and $N_\tau$. Usually the parameter combination is succinctly incorporated 
into the parameter name as, e.g. `l328f21b6390m00257m0694`. 

The `latticeParams` class collects all this information into one object, and contains some useful 
methods, for example the ability to convert one of these parameter strings into its corresponding 
float value. To instantiate this object, use
```Python
lp = latticeParams(Ns, Nt, beta, mass_l='', mass_s='')
```
Here `beta` and the quark masses are input strings coming from `l328f21b6390m00257m0694`. You can 
then, for instance get the float value corresponding to the string `6390` using
```Python
lp.getBeta()
```
There is also a method to get the `l328f21b6390m00257m0694` string, `getcparams`. You can also get 
the lattice spacing in [fm] and temperature in [MeV] using `geta` and `getT`, respectively. 
Finally `paramSummary()` prints a nice summary of all parameters going into your calculation 
to the screen.

There are always quark mass tables, which are implemented as dictionaries. Generally when we do 
our lattice calculations, we move along a line of constant physics where the ratio $m_s/m_l$ is 
fixed, and $m_s$ is fixed at its physical value. Hence if you know this ratio, $\beta$, and $N_\tau$, 
both quark masses are determined. For instance a small $N_\tau=12$ dictionary at the time of 
this writing is
```Python
class quarkMassNt12:

  """Lookup tables for Nt=12, Nf=2+1 flavors."""

  Table27={ '6794': ['00167','0450'],
            '6850': ['00157','0424'],
            '6910': ['00148','0401'] }
```
