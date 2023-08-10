# Statistical Physics

This module collects some basic methods and classes related to statistical physics.

## Critical exponents

There is a base `UniversalityClass` object from which all universality classes inherit. It holds critical exponents
as attributes along with an `exponentSummary()` method that prints to screen all critical exponents. If you want to
access 3-$d$, $\mathbb{Z}_2$ critical exponents, you can simply use, e.g.
```Python
univ = Z2_3d
univ.alpha
```

At present one finds the universality classes
- 3-$d$, O$(2)$
- 3-$d$, O$(3)$
- 3-$d$, O$(4)$
- 4-$d$, $\mathbb{Z}_2$
- 3-$d$, $\mathbb{Z}_2$
- 2-$d$, $\mathbb{Z}_2$

