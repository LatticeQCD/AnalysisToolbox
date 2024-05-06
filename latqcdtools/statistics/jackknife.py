#
# jackknife.py
#
# H. Dick, D. Clarke
#


import numpy as np
from latqcdtools.statistics.statistics import std_mean, std_err
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import parallel_function_eval
from latqcdtools.base.check import checkType
from latqcdtools.base.utilities import isHigherDimensional


def _pseudobins(jackknifeBins,avg):
    """ Calculate the 'pseudovalue' from the ith jackknife estimator. The pseudovalue is unbiased
    up to O(1/N), where N is the number of data. See e.g. eq. (1.1) of Miller, Biometrika 1974. """ 
    nblocks=len(jackknifeBins)
    return nblocks*np.array(avg) - (nblocks-1)*jackknifeBins


def _pareAxis(data,axis,nblocks):
    """ In order to ensure that all jackknife blocks have the same length, we pare the data along
    the conf_axis, so that nblocks divides the length of data along conf_axis."""
    data=np.array(data)
    length = data.shape[axis]
    new_length = length - (length % nblocks)
    indices = [slice(None)] * data.ndim
    indices[axis] = slice(new_length)
    return data[tuple(indices)]


def jackknife(f, data, numb_blocks=20, conf_axis=1, nproc=1, return_sample=False, args=()):
    """ Carry out a jackknife of an arbitrary function f of some data.

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
    """
    checkType(numb_blocks,int)
    checkType(conf_axis,int)
    checkType(nproc,int)
    checkType(return_sample,bool)
    if type(data) is tuple:
        logger.TBError('Does not support tuple data. Use latqcdtools.legacy.jackknife if you need this.')
    if numb_blocks <= 1:
        logger.TBError('Need numb_blocks > 1. Set numb_blocks=len(data) for remove-1 jackknife.')
    if not isHigherDimensional(data):
        conf_axis=0
    data = _pareAxis(data,conf_axis,numb_blocks)
    n = data.shape[conf_axis]
    total = f(data, *args)
    block_id = np.linspace(0, numb_blocks, n, endpoint=False).astype(np.int32)
    def fJ(i):
        sub_data = np.compress((block_id != i), data, axis=conf_axis)
        return f(sub_data, *args)
    J = np.array(parallel_function_eval(fJ,range(numb_blocks),nproc=nproc))
    bins = _pseudobins(J,total)
    if return_sample:
        return bins, std_mean(bins), std_err(bins) 
    else:
        return std_mean(bins), std_err(bins) 
