latqcdtools.base.initialize
=============

```python
TBRNG(seed=None) -> numpy.random._generator.Generator:
"""
The toolbox random number generator. This is O'Neill's PCG algorithm, which at the time of writing
is what numpy's default_rng() uses. We pin PCG64 explicitly, so the random numbers stay the same
even if numpy changes its default.

Args:
    seed (int or list of int, optional): Seed. Defaults to None, which chooses a random seed.

Returns:
    np.random.Generator
"""
```
```python
finalize():
"""
Some common tasks to do when you're done. 
"""
```
```python
gitHash() -> str:
"""
Obtain the current git hash. This assumes the Toolbox has been correctly installed.

Returns:
    str: git hash 
"""
```
```python
initialize(logFile=None):
"""
Some common tasks to do at the start of a run where you want to keep track of things. 
"""
```
```python
introduceYourself():
"""
Corporate branding. ASCII generated from https://patorjk.com. 
"""
```
