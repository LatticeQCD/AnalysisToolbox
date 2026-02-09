latqcdtools.physics.diracFreespectra
=============

```Python
_momenta_1d(L: int, antiperiodic: bool = False) -> numpy.ndarray:
'''
Return exactly L momentum modes.
n = -L//2, ..., L//2-1
p = 2π (n + shift)/L, shift=1/2 for antiperiodic
Works for even/odd L.
'''
```
```Python
class DiracOp(Lx, Ly, Lz, Lt, fermion='Wilson', bc_t='anti'):
```
```Python
class GammaMatrix():
'''The 4x4 gamma matrices used in Euclidean QFT.'''
```
