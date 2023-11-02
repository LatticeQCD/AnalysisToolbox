#
# jackknife.py
#
# H. Sandmeyer, D. Clarke
#
# A parallelized jackknife routine that can handle arbitrary return values of functions.
#

import numpy as np
import math
from latqcdtools.statistics.statistics import std_mean, std_err, meanArgWrapper
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import DEFAULTTHREADS, parallel_function_eval
from latqcdtools.base.utilities import isHigherDimensional
from latqcdtools.base.check import checkType


def _pseudo(mean, mean_i, numb_blocks):
    """ Calculate the 'pseudovalue' from the ith jackknife estimator. The pseudovalue is unbiased
    up to O(1/N), where N is the number of data. See e.g. eq. (1.1) of Miller, Biometrika 1974. """ 
    return numb_blocks*mean - (numb_blocks - 1)*mean_i


def _pseudo_val(mean, mean_i, numb_blocks):
    if type(mean) is tuple:
        retvalue = ()
        for i in range(len(mean)):
            retvalue += (_pseudo(np.array(mean[i]), np.array(mean_i[i]), numb_blocks),)
        return retvalue
    else:
        return _pseudo(np.array(mean), np.array(mean_i), numb_blocks)


def _getBlocksize(ndat,nblocks):
    if ndat < nblocks:
        logger.TBError('More blocks than data! nblocks, ndat =',nblocks,ndat)
    elif ndat==nblocks:
        return 1
    else:
        return math.floor(ndat/nblocks)


class nimbleJack:
    """ Class allowing for parallelization of the jackknife function. """

    def __init__(self, func, data, nblocks, confAxis, return_sample, args, nproc):

        checkType(nblocks,int)
        checkType(confAxis,int)
        checkType(return_sample,bool)
        checkType(nproc,int)
        self._func = func
        self._data = np.array(data)
        if nblocks==1:
            self._nblocks = len(data) 
        else:
            self._nblocks = nblocks
        if not isHigherDimensional(self._data):
            self._confAxis = 0
        else:
            self._confAxis = confAxis
        self._return_sample = return_sample
        self._args = args
        self._nproc = nproc

        # It could be that different observables have different numbers of measurements/configurations. We only
        # allow for that possibility by setting confAxis==1. If your configurations are ragged and along a
        # different axis, you're going to have to handle that yourself for now.
        if self._confAxis == 1:
            self._lengths = [len(self._data[i]) for i in range(len(self._data))]
            self._blocksizes = [_getBlocksize(length, self._nblocks) for length in self._lengths]
        else:
            self._length = self._data.shape[self._confAxis]
            self._blocksize = _getBlocksize(self._length,self._nblocks)

        if self._confAxis == 1:
            self._numb_observe = len(self._data)
            used_data = [self._data[i][:self._nblocks * self._blocksizes[i]] for i in range(self._numb_observe)]
            used_data = np.array(used_data)
        else:
            if self._length < self._nblocks:
                logger.TBError("More blocks than data. length, nblocks=", self._length, self._nblocks)
            swap_data = np.swapaxes(self._data, axis1=self._confAxis, axis2=0)
            used_data = np.swapaxes(swap_data[:self._blocksize*self._nblocks], 0, self._confAxis)
            if self._nblocks==self._length:
                if not np.array_equal(used_data,self._data):
                    logger.TBError('Something went wrong constucting used_data for remove-1 jackknife.')

        # Get initial estimator of the mean
        self._mean = meanArgWrapper(func, used_data, args)
        self._blockval = parallel_function_eval(self.getJackknifeEstimator, range(self._nblocks), nproc=self._nproc)
        self._mean = std_mean(self._blockval)
        self._error = std_err(self._blockval)

    def __repr__(self) -> str:
            return "jackknife"
        
    def getJackknifeEstimator(self, i):
        """ Gets the ith jackknife estimator from throwing away jackknife block i. """
        block_data = []
        if self._confAxis == 1:
            for k in range(self._numb_observe):
                block_data.append(self._data[k][0:i * self._blocksizes[k]])
                block_data[k] = np.concatenate((block_data[k], self._data[k][(i + 1) * self._blocksizes[k]:self._nblocks * self._blocksizes[k]]))
        else:
            # The idea is that the Markov time measurements are along the _confAxis. We construct the block_data
            # by swapping the conf_data to the 0th axis, so that we can cut out the ith configuration block
            # by just appending together using syntax like "swap_data[j]". Then we swap back at the end to restore
            # the original data structure.
            swap_data = np.swapaxes(self._data,axis1=self._confAxis,axis2=0)
            block_data = []
            for j in range(i*self._blocksize):
                block_data.append(swap_data[j])
            for j in range((i+1)*self._blocksize, self._nblocks*self._blocksize):
                block_data.append(swap_data[j])
            block_data = np.swapaxes(np.array(block_data), axis1=0, axis2=self._confAxis)
        mean_i = meanArgWrapper(self._func, block_data, self._args)
        return _pseudo_val(self._mean, mean_i, self._nblocks)

    def getResults(self):
        if self._return_sample:
            return self._blockval, self._mean, self._error
        else:
            return self._mean, self._error


def jackknife(func, data, numb_blocks=20, conf_axis=1, return_sample=False, args=(), nproc=DEFAULTTHREADS):
    """Jackknife routine for arbitray functions. This routine creates the jackknife like blocked subsets of data and
    passes them to the function in the same format as in the input data. So the idea is to write a function that
    computes an observable from a given data set. This function can be put into this jackkife routine and will get the
    jackknifed blocked data as input. Based on the output of the function, the jackkife mean and error are computed.
    The function may return multiple observables that are either scalars or numpy objects. You can pass a
    multidimensional object as data, but the jackknife function has to know from which axis blocks should be removed.
    This is controlled by conf_axis (default = 0 for one dimensional arrays and default = 1 for higher order arrays).
    Look into the __init__ function of the nimbleJack class to see in detail how this is implemented.

        Parameters
        ----------
        func : callable
            The function that calculates the observable

        data : array_like
            Input data

        numb_blocks : integer
            Number of jackknife blocks

        conf_axis : integer, optional, default = 0 for dim(data) = 1 and default = 1 for dim(data) >= 2
            Axis that should be resampled.

        return_sample : boolean, optional, default = False                                           
            Along with the mean and the error also return the results from the individual samples

        args : array_like or dict, default = ()
            optional arguments to be passed to func. If a dictionary they are passed as **args.

        nproc : integer
            Number of threads to use if you choose to parallelize. nproc=1 turns off parallelization.
    """
    jk = nimbleJack(func, data, numb_blocks, conf_axis, return_sample, args, nproc)
    return jk.getResults()


