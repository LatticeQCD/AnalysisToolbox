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
lp = latticeParams(Ns, Nt, beta, mass1='00257', mass2='0694', scaleType='fk', paramYear=2021, Nf='21')
```
Here `beta` and the quark masses are input strings coming from `l328f21b6390m00257m0694`. You can 
then, for instance, access the float value corresponding to the string `6390` via
```Python
lp.beta
```
The `scaleType` argument lets you choose what reference scale to use (`'fk'`, `'r1'`, or `'r0'`). How these
scales vary in lattice units as a function of $\beta$ has been calculated by the HotQCD collaboration; the
`paramYear` option lets you pick the year this function was computed, and `scaleYear` selects the year of the
physical-units scale determination.

Scale setting separates three distinct "worlds":
- `Nf`: the sea-quark content of your ensemble (`'21'`, `'111'`, `'211'`, `'2'`, `'3'`, `'5'`, or `None` for
  pure gauge). For `scaleType='r0'`, `Nf='21'` defaults `paramYear` to the `2012` $N_f=2+1$ parameterization,
  while `Nf=None` defaults to the `2017` pure-gauge parameterization.
- `aWorld`: the world in which the chosen $a(\beta)$ parameterization was fitted (e.g., `'Nf21'` or `'Nf0'`),
  resolved automatically via `get_aWorld(scaleType, year)`. If `aWorld` does not match your ensemble's `Nf`,
  a warning is logged (which can be suppressed with `ignoreWorldWarning()`).
- `physWorld`: the world from which the physical scale value in fm or MeV is taken. For `r0` and `r1`, only the
  `'Nf21'` determination is available. For `scaleType='fk'`, you can select `physWorld` among `'nature'` (PDG),
  `'Nf2'`, `'Nf21'`, and `'Nf211'` (defaulting to `'nature'` when `scaleYear` is `2012` or `2018`, and `'Nf21'`
  otherwise).

There is also a method to get the `l328f21b6390m00257m0694` string, `getcparams`. You can also get 
the lattice spacing in [fm] and temperature in [MeV] using `geta` and `getT`, respectively. 
Finally `paramSummary()` prints a nice summary of all parameters—including `Nf`, `aWorld`, and `physWorld`—going
into your calculation to the screen.

