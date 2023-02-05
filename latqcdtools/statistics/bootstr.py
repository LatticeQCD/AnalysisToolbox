#
# bootstr.py
#
# H. Sandmeyer, H.-T. Shu
#
# A parallelized bootstrap routine that can handle arbitrary return values of functions.
#
import numpy as np
from numpy.random import randint, normal, multivariate_normal
from latqcdtools.statistics.statistics import std_mean, std_dev, std_median, dev_by_dist
import concurrent.futures


# Recursive function to fill the sample
def recurs_append(data, sample_data, axis, conf_axis, sample_size, same_rand_for_obs, i, my_seed):

    np.random.seed(my_seed+i)
    if axis + 1 == conf_axis:
        numb_observe = len(data)
        if sample_size == 0:
            sample_sizes = [ len(i) for i in data ]
        else:
            sample_sizes = [sample_size]*len(data)

        if not same_rand_for_obs:
            randints = [randint(0, len(data[x]), size=sample_sizes[x]) for x in range(numb_observe)]
        else:
            tmp_rand = randint(0, len(data[0]), size=sample_sizes[0])
            randints = [tmp_rand]*numb_observe
        for x in range(numb_observe):
            sample_data.append(np.array(data[x])[randints[x]])
        return

    else:
        for i in range(len(data)):
            sample_data.append([])
            recurs_append(data[i], sample_data[i], axis + 1, conf_axis, sample_size, same_rand_for_obs, my_seed)


class nimbleBoot:

    def __init__(self, func, data, numb_samples, sample_size, same_rand_for_obs, conf_axis, return_sample, seed,
                 err_by_dist, args, parallelize, nproc):

        self._func=func
        self._data=np.array(data)
        self._numb_samples=numb_samples
        self._sample_size=sample_size
        self._same_rand_for_obs=same_rand_for_obs
        self._conf_axis=conf_axis
        self._return_sample=return_sample
        self._seed=seed
        self._err_by_dist=err_by_dist
        self._args=args

        if  self._data.ndim == 1:
            self._conf_axis = 0

        if self._seed is not None:
            my_seed = self._seed
        else:
            my_seed = randint(0, 7271978) #generate a time dependent seed
        if parallelize:
            sampleList=range(self._numb_samples)
            with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
                sampleval_i=executor.map(self.getBootstrapEstimator, sampleList,[my_seed]*self._numb_samples)
            self._sampleval=list(sampleval_i)
        else:
            sampleval=[]
            for i in range(self._numb_samples):
                sampleval.append(self.getBootstrapEstimator(i,my_seed))
            self._sampleval=sampleval

        if not self._err_by_dist:
            self._mean = std_mean(self._sampleval)
            self._error = std_dev(self._sampleval)
        else:
            self._mean = std_median(self._sampleval)
            self._error = dev_by_dist(self._sampleval)

    def getBootstrapEstimator(self,i,my_seed):
        sample_data = []
        if self._conf_axis == 0: # Case of one dimensional array is special
            np.random.seed(my_seed+i)
            if self._sample_size == 0:
                sample_size_tmp = len(self._data)
            else:
                sample_size_tmp = self._sample_size
            randints = randint(0, len(self._data), size=sample_size_tmp)
            sample_data = self._data[randints]

        else:
            axis = 0
            recurs_append(self._data, sample_data, axis, self._conf_axis, self._sample_size, self._same_rand_for_obs, i, my_seed)

        sample_data = np.array(sample_data)

        if isinstance(self._args, dict):
            mean_i = self._func(sample_data, **self._args)
        else:
            mean_i = self._func(sample_data, *self._args)
        return mean_i

    def getResults(self):
        if self._return_sample:
            return self._sampleval, self._mean, self._error
        else:
            return self._mean, self._error


def bootstr(func, data, numb_samples, sample_size = 0, same_rand_for_obs = False, conf_axis = 1, return_sample = False,
            seed = None, err_by_dist = False, args=(), parallelize=True, nproc=32):
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

        parallelize : boolean, optional, default = True
            Should you make it parallel? You may need to turn it off for some reason.

        nproc : integer, optional, default = 32
            Number of threads to use if you choose to parallelize.
    """
    bts = nimbleBoot(func, data, numb_samples, sample_size, same_rand_for_obs, conf_axis, return_sample, seed,
                     err_by_dist, args, parallelize, nproc)
    return bts.getResults()


class nimbleGaussianBoot:

    def __init__(self, func, data, data_std_dev, numb_samples, sample_size, same_rand_for_obs, return_sample, seed,
                 err_by_dist, useCovariance, Covariance, args, parallelize, nproc, asym_err):

        self._func=func
        self._data=np.array(data)
        self._data_std_dev=np.array(data_std_dev)
        self._numb_samples=numb_samples
        self._sample_size=sample_size
        self._same_rand_for_obs=same_rand_for_obs
        self._return_sample=return_sample
        self._seed=seed
        self._err_by_dist=err_by_dist
        self._useCovariance=useCovariance
        self._Covariance=Covariance
        self._args=args

        self._numb_observe = len(data)

        if self._seed is not None:
            my_seed = self._seed
        else:
            my_seed = randint(0, 7271978) # generate a time dependent seed

        if parallelize:
            sampleList=range(self._numb_samples)
            with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
                sampleval_i=executor.map(self.getGaussianBootstrapEstimator, sampleList,[my_seed]*self._numb_samples)
            self._sampleval=list(sampleval_i)
        else:
            sampleval=[]
            for i in range(self._numb_samples):
                sampleval.append(self.getGaussianBootstrapEstimator(i,my_seed))
            self._sampleval=sampleval

        if not self._err_by_dist:
            self._mean = std_mean(self._sampleval)
            self._error = std_dev(self._sampleval)
        else:
            self._mean = std_median(self._sampleval)
            self._error = dev_by_dist(self._sampleval, return_both_q=asym_err)

    def getGaussianBootstrapEstimator(self,i,my_seed):
        sample_data = []
        if not self._same_rand_for_obs:
            for k in range(self._numb_observe):
                np.random.seed(my_seed+i+(k+1)*self._numb_samples)
                if not self._useCovariance:
                    if self._sample_size == 1:
                        sample_data.append(normal(self._data[k], self._data_std_dev[k])) # k is the Nt data
                    else:
                        sample_data.append(normal(self._data[k], self._data_std_dev[k], self._sample_size))
                else:
                    if self._sample_size == 1:
                        if self._Covariance is None:
                            sample_data.append(multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2))) # k is the Nt data
                        else:
                            sample_data.append(multivariate_normal(self._data[k], self._Covariance[k])) # k is the Nt data
                    else:
                        if self._Covariance is None:
                            sample_data.append(multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2), self._sample_size))
                        else:
                            sample_data.append(multivariate_normal(self._data[k], self._Covariance[k], self._sample_size))
        else:
            np.random.seed(my_seed+i)
            for k in range(self._numb_observe):
                if not self._useCovariance:
                    if self._sample_size == 1:
                        sample_data.append(normal(self._data[k], self._data_std_dev[k])) # k is the Nt data
                    else:
                        sample_data.append(normal(self._data[k], self._data_std_dev[k], self._sample_size))
                else:
                    if self._sample_size == 1:
                        if self._Covariance is None:
                            sample_data.append(multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2))) # k is the Nt data
                        else:
                            sample_data.append(multivariate_normal(self._data[k], self._Covariance[k])) # k is the Nt data
                    else:
                        if self._Covariance is None:
                            sample_data.append(multivariate_normal(self._data[k], np.diag(self._data_std_dev[k]**2), self._sample_size))
                        else:
                            sample_data.append(multivariate_normal(self._data[k], self._Covariance[k], self._sample_size))

        sample_data = np.array(sample_data)

        if isinstance(self._args, dict):
            mean_i = self._func(sample_data, **self._args)
        else:
            mean_i = self._func(sample_data, *self._args)

        return mean_i

    def getResults(self):
        if self._return_sample:
            return self._sampleval, self._mean, self._error
        else:
            return self._mean, self._error


def bootstr_from_gauss(func, data, data_std_dev, numb_samples, sample_size = 1, same_rand_for_obs = False,
                       return_sample = False, seed = None, err_by_dist = False, useCovariance = False,
                       Covariance = None, args = (), parallelize = True, nproc = 32, asym_err=False):
    """Same as standard bootstrap routine, but the data are generated by gaussian noise around the mean values in data.
    The width of the distribution is controlled by data_std_dev. Note, that the function has to average over samples.
    This means that data_std_dev should always be the standard deviation of a single measurement and not the standard
    deviation of a mean."""
    data = np.asarray(data)
    data_std_dev = np.asarray(data_std_dev)

    bts_gauss = nimbleGaussianBoot(func, data, data_std_dev, numb_samples, sample_size, same_rand_for_obs,
                                   return_sample, seed, err_by_dist, useCovariance, Covariance, args, parallelize,
                                   nproc, asym_err)
    return bts_gauss.getResults()
