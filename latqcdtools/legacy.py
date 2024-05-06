# 
# legacy.py                                                               
# 
# D. Clarke, H. Sandmeyer 
# 
# Here we collect "legacy code". This is code that has to meet some requirements:
#   1. It is hard for David to read,
#   2. and therefore hard for him to maintain,
#   3. while not really being used. 
# 
# This module exists for users who still want to utilize some of this legacy
# functionality. Before making something legacy code, please put some kind of
# warning or error somewhere.
# 


import numpy as np
import math
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.base.utilities import isHigherDimensional
from latqcdtools.base.speedify import parallel_function_eval, DEFAULTTHREADS
from latqcdtools.statistics.statistics import meanArgWrapper
from latqcdtools.base.initialize import DEFAULTSEED


# Suppose the wrapped function is `sum` and `data` is `((1, 2, 3), (4, 5, 6), (7, 8, 9))`. 
# Then reduce_tuple will pass `(1, 4, 7)`, `(2, 5, 8)`, `(3, 6, 9)` to the `sum` function 
# respectively and return `(sum(1, 4, 7), sum(2, 5, 8), sum(3, 6, 9))`.
def reduce_tuple(func):
    def func_wrapper(data, *args, **kwargs):
        if type(data[0]) is tuple:
            retvalue = ()
            for i in range(len(data[0])):
                obj_array = []
                for k in range(len(data)):
                    obj_array.append(data[k][i])
                retvalue += (func(obj_array, *args, **kwargs),)
            return retvalue
        else:
            return func(data, *args, **kwargs)
    return func_wrapper


@reduce_tuple
def std_median(data, axis = 0):
    return np.median(data, axis)


@reduce_tuple
def std_mean(data, axis = 0):
    return np.mean(data, axis)


@reduce_tuple
def std_var(data, axis = 0):
    data = np.asarray(data)
    return np.var(data, axis = axis, ddof = 1)


@reduce_tuple
def std_dev(data, axis = 0):
    data = np.asarray(data)
    return np.std(data, axis = axis, ddof = 1)


@reduce_tuple
def std_err(data, axis = 0):
    data = np.asarray(data)
    return std_dev(data, axis) / np.sqrt(data.shape[axis])


def _pseudo(mean, mean_i, numb_blocks):
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


@reduce_tuple
def dev_by_dist(data, axis=0, return_both_q=False, percentile=68):
    data = np.asarray(data)
    median = np.nanmedian(data, axis)
    numb_data = data.shape[axis]
    idx_dn = max(int(np.floor((numb_data-1) / 2 - percentile/100 * (numb_data-1) / 2)), 0)
    idx_up = min(int(np.ceil((numb_data-1) / 2 + percentile/100 * (numb_data-1) / 2)), numb_data-1)
    sorted_data = np.sort(data - np.expand_dims(median, axis), axis=axis)
    q_l = np.take(sorted_data, idx_dn, axis)
    q_r = np.take(sorted_data, idx_up, axis)
    if return_both_q:
        return np.abs(q_l), np.abs(q_r)
    else:
        return np.max(np.stack((np.abs(q_l), np.abs(q_r)), axis=0), axis=0)


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


def _autoSeed(seed):
    """ We use seed=None to flag the seed should be automatically chosen. The problem is that we need
    seed to be an integer when enforcing that different bootstrap samples use different seeds. """
    if seed is None:
        return np.random.randint(0,DEFAULTSEED)
    else:
        return seed


def recurs_append(data, sample_data, axis, conf_axis, sample_size, same_rand_for_obs, i, my_seed):
    """ Recursive function to fill the sample. """

    rng = np.random.default_rng(my_seed+i)
    if axis + 1 == conf_axis:
        numb_observe = len(data)
        if sample_size == 0:
            sample_sizes = [ len(j) for j in data ]
        else:
            sample_sizes = [sample_size]*len(data)

        if not same_rand_for_obs:
            randints = [rng.integers(0, len(data[x]), size=sample_sizes[x]) for x in range(numb_observe)]
        else:
            tmp_rand = rng.integers(0, len(data[0]), size=sample_sizes[0])
            randints = [tmp_rand]*numb_observe
        for x in range(numb_observe):
            sample_data.append(np.array(data[x])[randints[x]])
        return

    else:
        for j in range(len(data)):
            sample_data.append([])
            recurs_append(data[j], sample_data[j], axis + 1, conf_axis, sample_size, same_rand_for_obs, i, my_seed)


class nimbleBoot:

    def __init__(self, func, data, numb_samples, sample_size, same_rand_for_obs, conf_axis, return_sample, seed,
                 err_by_dist, args, nproc):

        self._func=func
        self._data=np.array(data)
        self._numb_samples=numb_samples
        self._sample_size=sample_size
        self._same_rand_for_obs=same_rand_for_obs
        self._conf_axis=conf_axis
        self._return_sample=return_sample
        self._seed=_autoSeed(seed)
        checkType(self._seed,int)
        self._err_by_dist=err_by_dist
        self._args=args
        self._nproc = nproc 

        if  self._data.ndim == 1:
            self._conf_axis = 0

        self._sampleval = parallel_function_eval(self.getBootstrapEstimator,range(self._numb_samples),nproc=self._nproc,args=(self._seed,))

        if not self._err_by_dist:
            self._mean = std_mean(self._sampleval)
            self._error = std_dev(self._sampleval)
        else:
            self._mean = std_median(self._sampleval)
            self._error = dev_by_dist(self._sampleval)

    def __repr__(self) -> str:
        return "nimbleBoot"

    def getBootstrapEstimator(self,i,my_seed):
        sample_data = []
        if self._conf_axis == 0: # Case of one dimensional array is special
            rng = np.random.default_rng(my_seed+i)
            if self._sample_size == 0:
                sample_size_tmp = len(self._data)
            else:
                sample_size_tmp = self._sample_size
            randints = rng.integers(0, len(self._data), size=sample_size_tmp)
            sample_data = self._data[randints]

        else:
            axis = 0
            recurs_append(self._data, sample_data, axis, self._conf_axis, self._sample_size, self._same_rand_for_obs, i, my_seed)

        sample_data = np.array(sample_data)

        return meanArgWrapper(self._func, sample_data, self._args)

    def getResults(self):
        if self._return_sample:
            return self._sampleval, self._mean, self._error
        else:
            return self._mean, self._error


def bootstr(func, data, numb_samples, sample_size = 0, same_rand_for_obs = False, conf_axis = 1, return_sample = False,
            seed = None, err_by_dist = False, args=(), nproc=DEFAULTTHREADS):
    """Bootstrap for arbitrary functions. This routine resamples the data and passes them to the function in the same
    format as in the input. So the idea is to write a function that computes an observable from a given data set. This
    function can be put into this bootstrap routine and will get bootstrap samples as input. Based on the output of the
    function, the bootstrap mean and error are computed. The function may return multiple observables that are either
    scalars or numpy objects. You can pass a multidimensional object as data, but the bootstrap function has to know
    which axis should be resampled which is controlled by conf_axis (default = 0 for one dimensional arrays and
    default = 1 for higher order arrays.)

        Parameters
        ----------
        func : callable
            The function that calculates the observable

        data : array_like
            Input data

        numb_samples : integer
            Number of bootstrap samples

        sample_size : integer, optional, default = 0
            Size of sample

        same_rand_for_obs : boolean, optional, default = False
            Use the same random numbers for each observable accessed by index conf_axis - 1. Please note:
                - Objects that are accessed by an axis >= conf_axis + 1 do always have the same random numbers.
                - Objects that are accessed by axis conf_axis < conf_axis - 1 never share the same random numbers.

        conf_axis : integer, optional, default = 0 for dim(data) = 1 or default = 1 for dim(data) >= 2
            Axis that should be resampled

        return_sample : boolean, optional, default = False
            Along with the mean and the error also return the results from the individual samples

        seed: integer, optional, default = None
            seed for the random generator. If None, the default seed from numpy is used (probably from time)

        same_rand_for_obs : boolean, optional, default = False
            same random numbers per observable

        err_by_dist : boolean, optional, default = False
            Compute the error from the distribution using the median and the 68% quantile

        args : array_like or dict, default = ()
            optional arguments to be passed to func. If a dictionary the are passed as **args.

        nproc : integer
            Number of threads to use if you choose to parallelize. nproc=1 turns off parallelization.
    """
    bts = nimbleBoot(func, data, numb_samples, sample_size, same_rand_for_obs, conf_axis, return_sample, seed,
                     err_by_dist, args, nproc)
    return bts.getResults()


class nimbleGaussianBoot:

    def __init__(self, func, data, data_std_dev, numb_samples, sample_size, same_rand_for_obs, return_sample, seed,
                 err_by_dist, useCovariance, Covariance, args, nproc, asym_err):

        checkType(numb_samples,int)
        checkType(sample_size,int)
        checkType(same_rand_for_obs,bool)
        checkType(return_sample,bool)
        checkType(err_by_dist,bool)
        checkType(useCovariance,bool)
        checkType(nproc,int)
        checkType(asym_err,bool)
        self._func=func
        self._data=np.array(data)
        self._data_std_dev=np.array(data_std_dev)
        self._numb_samples=numb_samples
        self._sample_size=sample_size
        self._same_rand_for_obs=same_rand_for_obs
        self._return_sample=return_sample
        self._seed=_autoSeed(seed)
        checkType(self._seed,int)
        self._err_by_dist=err_by_dist
        self._useCovariance=useCovariance
        self._Covariance=Covariance
        self._args=args
        self._numb_observe = len(data)
        self._nproc = nproc 

        self._sampleval=parallel_function_eval(self.getGaussianBootstrapEstimator,range(self._numb_samples),nproc=self._nproc,args=(self._seed,))

        if not self._err_by_dist:
            self._mean = std_mean(self._sampleval)
            self._error = std_dev(self._sampleval)
        else:
            self._mean = std_median(self._sampleval)
            self._error = dev_by_dist(self._sampleval, return_both_q=asym_err)

    def __repr__(self) -> str:
        return "nimbleGaussianBoot"

    def getGaussianBootstrapEstimator(self,i,my_seed):
        sample_data = []
        if not self._same_rand_for_obs:
            for k in range(self._numb_observe):
                rng=np.random.default_rng(my_seed+i+(k+1)*self._numb_samples)
                if not self._useCovariance:
                    if self._sample_size == 1:
                        sample_data.append(rng.normal(self._data[k], self._data_std_dev[k]))
                    else:
                        sample_data.append(rng.normal(self._data[k], self._data_std_dev[k], self._sample_size))
                else:
                    if self._sample_size == 1:
                        if self._Covariance is None:
                            sample_data.append(rng.multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2)))
                        else:
                            sample_data.append(rng.multivariate_normal(self._data[k], self._Covariance[k]))
                    else:
                        if self._Covariance is None:
                            sample_data.append(rng.multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2), self._sample_size))
                        else:
                            sample_data.append(rng.multivariate_normal(self._data[k], self._Covariance[k], self._sample_size))
        else:
            rng=np.random.default_rng(my_seed+i)
            for k in range(self._numb_observe):
                if not self._useCovariance:
                    if self._sample_size == 1:
                        sample_data.append(rng.normal(self._data[k], self._data_std_dev[k]))
                    else:
                        sample_data.append(rng.normal(self._data[k], self._data_std_dev[k], self._sample_size))
                else:
                    if self._sample_size == 1:
                        if self._Covariance is None:
                            sample_data.append(rng.multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2)))
                        else:
                            sample_data.append(rng.multivariate_normal(self._data[k], self._Covariance[k]))
                    else:
                        if self._Covariance is None:
                            sample_data.append(rng.multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2), self._sample_size))
                        else:
                            sample_data.append(rng.multivariate_normal(self._data[k], self._Covariance[k], self._sample_size))

        sample_data = np.array(sample_data)

        return meanArgWrapper(self._func,sample_data,self._args) 

    def getResults(self):
        if self._return_sample:
            return self._sampleval, self._mean, self._error
        else:
            return self._mean, self._error


def bootstr_from_gauss(func, data, data_std_dev, numb_samples, sample_size = 1, same_rand_for_obs = False,
                       return_sample = False, seed = None, err_by_dist = True, useCovariance = False,
                       Covariance = None, args = (), nproc = DEFAULTTHREADS, asym_err=False):
    """Same as standard bootstrap routine, but the data are generated by gaussian noise around the mean values in data.
    The width of the distribution is controlled by data_std_dev. Note, that the function has to average over samples.
    This means that data_std_dev should always be the standard deviation of a single measurement and not the standard
    deviation of a mean."""
    data = np.asarray(data)
    data_std_dev = np.asarray(data_std_dev)

    bts_gauss = nimbleGaussianBoot(func=func, data=data, data_std_dev=data_std_dev, numb_samples=numb_samples, 
                                   sample_size=sample_size, same_rand_for_obs=same_rand_for_obs,return_sample=return_sample, 
                                   seed=seed, err_by_dist=err_by_dist, useCovariance=useCovariance, Covariance=Covariance, 
                                   args=args, nproc=nproc, asym_err=asym_err)
    return bts_gauss.getResults()
