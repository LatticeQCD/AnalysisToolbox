# Writing Tests

With each new module you write, please add a test. All the tests for the LatticeToolbox can be found in the
`testing` folder. Once you have written your test, add a line for it to `testing/runTests.bash`.
If you make any changes to the code, call
```shell
bash runTests.bash
```
from the `testing` folder to make sure nothing is broken.

There are a couple modules that assist with testing. The module
```Python
import latqcdtools.math.check
```
contains a method 
```Python
rel_check(a, b, prec=1e-6)
``` 
to compare the floats `a` and `b` up to precision `prec`. This is essentially
a wrapper for `math.isclose`, and it's what we recommend especially when `a=0` or `b=0`. 

One also finds in this module 
```Python
print_results(res, res_true, res_err = None, res_err_true = None, text = "", prec = 1e-4)
```
which is a convenient method for comparing two array-like objects `res` and `res_true`.

It is also sometimes important during testing to time your code. This can be important for example to verify that
parallelizing indeed speeds up your code. (Maybe something could go wrong here for an unexpected reason.) the module
```Python
import latqcdtools.base.utilities
```
includes the `timer` class. A call to
```Python
timer t
t.printTiming()
```
prints the time in seconds that elapsed since either the timer was instantiated or
since the last time `printTiming()` was called.
