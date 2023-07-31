# Speedify

The module
```Python
import latqcdtools.base.speedify
```
includes methods for easily making your Python scripts faster. In particular we would like
to draw your attention to the `@compile` decorator and the `parallel_function_eval` methods.
The former is a decorator that effectively wraps `@njit` from [Numba](https://numba.pydata.org/). 
The latter calls a parallelizer behind the scenes to accelerate a for-loop. At the moment we have
implemented only `concurrent.futures` and `pathos.pools`, but it is straightforward to
extend this to your favorite parallelizer.

Numba's just-in-time compilation takes a generally nontrivial amount of time, so whether
one compiles a method may depend on the context, e.g. how many times that method is called.
Hence we have implemented the `@compile` decorator. All methods in a module using `@compile`
can be flagged to compile or not using `numbaON()` and `numbaOFF()`. This allows you to
do targeted compilation in your code; for instance
```Python
from latqcdtools.base.speedify import numbaON, numbaOFF
numbaON()
# gaugeField contains some methods that can be compiled. If the lattice is small enough,
# this compilation may not be worth your time.
from latqcdtools.physics.gauge import gaugeField
numbaOFF()
```
Incidentally, `latqcdtools.physics.gauge` includes some examples of the `@compile` decorator
in action, in case you want to try it yourself.

Next we show an example of replacing a for-loop with `parallel_function_eval`. Consider
a function `myFunc` **that is itself not already parallelized**. Then
```Python
result = []
for item in myList:
    result.append(myFunc(item))
```
can be replaced with
```Python
result = parallel_function_eval(myFunc,myList)
```
The function prototype looks like
```Python
parallel_function_eval(function, input_array, args=(), nproc=DEFAULTTHREADS, parallelizer=DEFAULTPARALLELIZER)
```
The `DEFAULTPARALLELIZER` is `pathos.pools` because it plays more nicely with classes, but you can also
use `concurrent.futures` if you prefer. `DEFAULTTHREADS` is computed behind the scenes to use the maximum
number of threads available to your system, minus 2 (so that we don't hog all your resources).