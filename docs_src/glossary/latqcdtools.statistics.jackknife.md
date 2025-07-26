latqcdtools.statistics.jackknife
=============

```Python
_pareAxis(data, axis, nblocks):
'''
In order to ensure that all jackknife blocks have the same length, we pare the data along
the conf_axis, so that nblocks divides the length of data along conf_axis.
'''
```
```Python
_pseudobins(jackknifeBins, avg):
'''
Calculate the 'pseudovalue' from the ith jackknife estimator. The pseudovalue is unbiased
up to O(1/N), where N is the number of data. See e.g. eq. (1.1) of Miller, Biometrika 1974. 
'''
```
```Python
jackknife(f, data, numb_blocks=20, conf_axis=1, nproc=1, return_sample=False, args=()):
'''
Carry out a jackknife of an arbitrary function f of some data.

Args:
    f (func)
    data (array-like)
    numb_blocks (int, optional): Number of jackknife blocks. Defaults to 20.
    conf_axis (int, optional): The axis that represents the configuration axis, i.e., measurements
      are along this axis. Defaults to 1 when data has dimension of at least 2; else defaults to 0.
    nproc (int, optional): Number of threads to use. Defaults to 1. Increasing nproc is likely to
      have a deleterious impact on performance when f is a fast function. When f is slower, for
      example if you are doing curve fits or integrating in the block, increasing nproc may
      improve performance. Feel free to benchmark it for your use-case. 
    return_sample (bool, optional): Return pseudovalues? Defaults to False.
    args (tuple, optional)

Returns:
    jackknife mean and error (optionally pseudovalues) 
'''
```
