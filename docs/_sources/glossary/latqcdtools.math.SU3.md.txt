latqcdtools.math.SU3
=============

```Python
class SU3(mat=None):
'''
A member of the Lie group SU(3). Implemented as a subclass of the np.ndarray class. This gives us access already
to all the nice features of np.ndarray and lets us leverage the speed of numpy.
    g.trace()
    g.det()
    g.dagger()
    g[i,j], which can be used to access and assign
    g + h
    g*h = g@h
    2*g
'''
```
