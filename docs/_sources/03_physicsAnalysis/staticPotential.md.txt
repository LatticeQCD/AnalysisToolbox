# The zero-temperature static potential

In this module are mehtods related to the zero-temperature static potential. This potential is sometimes used to
renormalize Polyakov-loop observables, to extract reference scales like $r_0$, or to improve short-distance
lattice spacings at tree-level.

## Tree-level improved distances
Short distances calculated directly from multiplying the separation by the lattice 
spacing are plagued by lattice artifacts. Therefore if one is interested in how an observable depends on $r$, they
should use tree-level improved distances for small $r$. This calculation depends on $N_s$. Given a maximum
squared distance to improve `r2max`, the method
```Python
impdist(Ns,r2max)
```
returns a list of improved distances.

## Zero temperature quark potential
The method
```Python
V_Teq0(r)
```
takes a distance $r$ in [fm] and returns $V$ in [MeV]. The potential comes from a three-parameter Levenberg-Marquardt 
fit of the data in Fig. 14 [here](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.094503) to the 
form $V(r)=a/r+b*r+c$, which is the Cornell parameterization.

There are also fit functions like `fitV_Teq0`, which can be used to extract fit parameters using a Cornell ansatz.
Other possibilities include a parameterization including one-loop and two-loop corrections to the Coulomb part.

