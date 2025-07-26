latqcdtools.physics.runningCoupling
=============

```Python
CA(Nc) -> float:
'''
Casimir operator of adjoint representation of SU(Nc).

Args:
    Nc (int): Number of colors.

Returns:
    float : CA 
'''
```
```Python
CF(Nc) -> float:
'''
Casimir operator of defining representation of SU(Nc).

Args:
    Nc (int): Number of colors.

Returns:
    float : CF 
'''
```
```Python
b0(Nf, Nc=3) -> float:
'''
Universal one-loop beta coefficient for SU(Nc), expansion in g.  

Args:
    Nf (int): Number of active quark flavors. 
    Nc (int): Number of colors.

Returns:
    float : b0 
'''
```
```Python
b1(Nf, Nc=3) -> float:
'''
Universal two-loop beta coefficient for SU(Nc), expansion in g. 

Args:
    Nf (int): Number of active quark flavors.
    Nc (int): Number of colors.

Returns:
    float : b1 
'''
```
```Python
b2_dimreg_MSbar(Nf, Nc=3) -> float:
'''
Three-loop beta coefficient for SU(Nc) using dimensional regularization in MS-bar scheme. 

Args:
    Nf (int): Number of active quark flavors.
    Nc (int): Number of colors.

Returns:
    float : b2_MS-bar (dim reg) 
'''
```
```Python
b3_dimreg_MSbar(Nf, Nc=3) -> float:
'''
Four-loop beta coefficient for SU(Nc) using dimensional regularization in MS-bar scheme. 

Args:
    Nf (int): Number of active quark flavors.
    Nc (int): Number of colors.

Returns:
    float : b3_MS-bar (dim reg) 
'''
```
```Python
beta_func(beta, Nf=3) -> float:
'''
QCD asymptotic scaling relation to two loops.

Args:
    beta (float-like)
    Nf (int, optional): Number of active quark flavors. Defaults to 3.

Returns:
    float: f_as(beta) 
'''
```
