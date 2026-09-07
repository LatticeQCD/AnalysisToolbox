latqcdtools.math.pade
=============

```python
class singlePointPade(c, p, q, x0=0.0):
"""
Single-point Padé approximant. Given the Taylor coefficients c of some function f about
a point x0, i.e.
    f(x) = sum_k c[k] (x-x0)**k,
build the rational function
    R(x) = P(x-x0) / Q(x-x0),   deg(P) = p,  deg(Q) = q,
whose own Taylor expansion about x0 agrees with f through order p+q. The denominator is
normalized so that Q(0) = 1. The object is callable and is backed by a Rational.

The numerator degree p and denominator degree q must both be given. They have to satisfy
p >= 0, q >= 0, and p + q <= len(c) - 1, since matching the series through order p+q needs
the coefficients c[0..p+q].

Example:
    from math import factorial
    c = [1/factorial(k) for k in range(7)]   # Taylor coeffs of exp about 0
    R = singlePointPade(c, p=3, q=3)
    R(0.5)                                    # ~ exp(0.5)
"""
```
