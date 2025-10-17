# SU(N) matrices 

In `latqcdtools.math.SU3` one finds the `SU3` class, which inherits from the `numpy.matrix` class. This way you
can get all the nice, fast features of `numpy`, like already defined matrix multiplication, transposes, and so on.

The `SU3` class also has some extra functionality, most importantly given the first two rows, the method
`su3unitarize` will use unitarity to compute the last one. This allows reconstruction of compressed configurations
in the `confReader` class, described [here](../interfacing/confLoader.md). The unitarization can be compiled
with `numbaON()`, as described [here](../base/speedify.md).

You can instantiate an `SU(3)` matrix object as
```Python
U = SU3()
```
This object is needed especially for `gaugeField` objects, described [here](../physicsAnalysis/gauge.md).
Some other implemented methods special for `SU3` objects include
- `setToRandom`: Sets to a random group element.
- `isSU3`: Checks that it has determinant 1 and is unitary.

Similarly there is some support for `SU(2)` objects in `latqcdtools.math.SU2`.
