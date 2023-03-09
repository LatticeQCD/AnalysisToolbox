# SU(3) matrices 

In `latqcdtools/math/SU3.py` one finds the `SU3` class, which inherits from the `numpy.matrix` class. This way you
can get all the nice, fast features of `numpy`, like already defined matrix multiplication, transposes, and so on.

The `SU(3)` class also has some extra functionality, most importantly given the first two rows, the method
`su3unitarize` will use unitarity to compute the last one. This allows reconstruction of compressed configurations
in the `confReader` class, described [here](../interfacing/confLoader.md).

You can instantiate an `SU(3)` matrix object as
```Python
U = SU3()
```
This object is needed especially for `gaugeField` objects, described [here](../physicsAnalysis/gauge.md).
