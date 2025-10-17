# Statistical Physics

This module collects some basic methods and classes related to statistical physics.
These can be imported from `latqcdtools.physics.statisticalPhysics`.

## Critical exponents

There is a base `UniversalityClass` object from which all universality classes inherit. It holds critical exponents
as attributes along with an `exponentSummary()` method that prints to screen all critical exponents. If you want to
access 3-$d$, $\mathbb{Z}_2$ critical exponents, you can simply use, e.g.
```Python
univ = Z2_3d
univ.alpha
```

## Reweighter

There is also a basic reweighter. Here we give a basic example to reweight a magnetic susceptibility:

```Python
def RWSUSC(data,xRW,x0) -> float:
    """ Reweight the susceptibility. The susceptibility is an observable that is
    defined in terms of expectation values. At the same time, we think of the
    reweight() method as a redefined expectation value.

    Args:
        data (list): a list [M, E] 
        xRW (float): the point we are RWing to 
        x0 (float): the starting point (plays role 1/T)

    Returns:
        float: reweighted susceptibility 
    """
    X = data[0]
    S = data[1]
    return x0*( reweight(X**2,xRW,x0,S) - reweight(X,xRW,x0,S)**2 )
```
