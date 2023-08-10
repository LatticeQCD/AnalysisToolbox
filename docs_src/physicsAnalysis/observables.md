# Polyakov loop observables


The Polyakov loop $P$ gives the order parameter of the deconfinement transition in the quenched limit. The imaginary
part is the order parameter of the Roberge-Weiss transition. It is also occasionally of interest at finite quark
mass and real $\mu$. You can find a bunch of Polyakov loop related observables in the `polyakovTools` module.

Initialize a `polyakovTools` object with
```Python
pt = polyakovTools(Ns, Nt)
```
with spatial extent `Ns` and Euclidean time extent `Nt`. One can then get $\langle |P|\rangle$ from
```Python
numpy.mean( pt.absP(RePmeas, ImPmeas) )
```
where `RePmeas` and `ImPmeas` are `numpy` arrays of measurements of real and imaginary parts of $P$, respectively. 