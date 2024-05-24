# Writing Tests

With each new module you write, please add a test. All the tests for the AnalysisToolbox can be found in the
`tests` folder.
If you make any changes to the code, call
```shell
pytest
```
from the top level. 

There are a couple modules that assist with testing. The module
```Python
import latqcdtools.math.math
```
contains a method 
```Python
rel_check(a, b, prec=1e-6)
``` 
to compare the floats `a` and `b` up to precision `prec`. This is essentially
a wrapper for `math.isclose`, and it's what we recommend especially when `a=0` or `b=0`. 

In the `latqcdtools.testing` module, one finds
```Python
print_results(res, res_true, res_err = None, res_err_true = None, text = "", prec = 1e-10, abs_prec=None)
```
which is a convenient method for comparing two array-like objects `res` and `res_true`, optionally allowing
you to compare errors as well. This returns `True` if `res` and `res_true` agree. Most test scripts are
structured to summarize the results of all tests using `concludeTest()`, which is also in this module,
and prints a standard message depending on the test outcomes. For instance
```Python
lpass = True
lpass *= print_results( res1, res_true1 )
lpass *= print_results( res2, res_true2 )
concludeTest(lpass)
```
Finally if you would like your test to be statistical, you can do Z-tests between two vectors given
their Gaussian error bars using
```Python
gaudif_results(res, res_err, res_true, res_err_true, text = "", qcut=0.05, testMode=True)
```
You can read more about our Z-test implementation [here](../dataAnalysis/statistics.md)


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
