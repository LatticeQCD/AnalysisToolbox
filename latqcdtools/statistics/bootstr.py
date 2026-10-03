#
# bootstr.py
#
# H. Sandmeyer, H.-T. Shu
#
# A parallelized bootstrap routine that can handle arbitrary return values of functions.
#


import numpy as np
from latqcdtools.statistics.statistics import meanArgWrapper, std_mean, std_dev, dev_by_dist
from latqcdtools.base.speedify import DEFAULTTHREADS, parallel_function_eval
from latqcdtools.base.initialize import DEFAULTSEED, TBRNG
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger


def _autoSeed(seed) -> int:
    """ 
    We use seed=None to flag the seed should be automatically chosen. The problem is that we need
    seed to be an integer when enforcing that different bootstrap samples use different seeds. 
    """
    if seed is None:
        return int(TBRNG().integers(0,DEFAULTSEED))
    else:
        checkType('int',seed=seed)
        return seed


def _biasCorrect(fbar,sampleval):
    """ 
    Bias-corrected central value 2 f(xbar) - <f*>. Here f(xbar) is func of the full data, <f*> is the
    average over bootstrap samples, and <f*> - f(xbar) estimates the O(1/N) bias of f(xbar). The result
    carries a Monte Carlo error of about (bootstrap error)/sqrt(numb_samples), even for linear func. 
    """
    return 2*np.asarray(fbar) - std_mean(sampleval)


class nimbleBoot:

    def __init__(self, func, data, numb_samples, sample_size, same_rand_for_obs, conf_axis, return_sample, seed,
                 err_by_dist, args, nproc):

        checkType('int',numb_samples=numb_samples)
        checkType(bool,same_rand_for_obs=same_rand_for_obs)
        checkType(bool,return_sample=return_sample)
        checkType(bool,err_by_dist=err_by_dist)
        checkType('int',nproc=nproc)
        self._func=func
        try:
            self._data=np.array(data)
        except ValueError:
            logger.TBRaise('All observables must have the same number of configurations.')
        self._numb_samples=numb_samples
        if sample_size is not None:
            checkType('int',sample_size=sample_size)
        self._sample_size=sample_size
        self._same_rand_for_obs=same_rand_for_obs
        self._conf_axis=conf_axis
        self._return_sample=return_sample
        self._seed=_autoSeed(seed)
        self._err_by_dist=err_by_dist
        self._args=args
        self._nproc = nproc 

        if self._data.ndim == 1:
            self._conf_axis = 0
        if self._conf_axis >= self._data.ndim:
            logger.TBRaise('conf_axis',self._conf_axis,'out of range for data with ndim',self._data.ndim)

        self._sampleval = parallel_function_eval(self.getBootstrapEstimator,range(self._numb_samples),nproc=self._nproc,args=(self._seed,))

        self._mean = _biasCorrect(meanArgWrapper(self._func, self._data, self._args), self._sampleval)
        if not self._err_by_dist:
            self._error = std_dev(self._sampleval)
        else:
            self._error = dev_by_dist(self._sampleval)

    def __repr__(self) -> str:
        return "nimbleBoot"

    def getBootstrapEstimator(self,i,my_seed):
        # Seeding with the pair (seed, i) gives unrelated streams for different seeds and samples.
        rng = TBRNG([my_seed,i])
        nconf = self._data.shape[self._conf_axis]
        if self._sample_size is None:
            sample_size = nconf
        else:
            sample_size = self._sample_size

        # Every index before conf_axis labels an observable.
        obs_shape = self._data.shape[:self._conf_axis]
        randints  = np.empty(obs_shape + (sample_size,), dtype=int)
        if self._same_rand_for_obs:
            # One draw of configurations, shared by all observables (perfectly correlated)
            randints[...] = rng.integers(0, nconf, size=sample_size)
        else:
            # An independent draw of configurations for each observable
            for obs in np.ndindex(obs_shape):
                randints[obs] = rng.integers(0, nconf, size=sample_size)

        # Indices after conf_axis belong to one configuration, so they share its random index.
        randints    = randints.reshape(randints.shape + (1,)*(self._data.ndim - self._conf_axis - 1))
        sample_data = np.take_along_axis(self._data, randints, axis=self._conf_axis)

        return meanArgWrapper(self._func, sample_data, self._args)

    def getResults(self):
        if self._return_sample:
            return np.array(self._sampleval), self._mean, self._error
        else:
            return self._mean, self._error


def bootstr(func, data, numb_samples, sample_size = None, same_rand_for_obs = False, conf_axis = 1, return_sample = False,
            seed = None, err_by_dist = True, args=(), nproc=DEFAULTTHREADS):
    """
    Bootstrap for arbitrary functions. This routine resamples the data and passes them to func in the same
    format as the input, so func should compute an observable from a given data set. The central value is
    the bias-corrected estimate 2 f(xbar) - <f*>, where f(xbar) is func of the full data and <f*> is the
    average over bootstrap samples. Like the jackknife, this removes the O(1/N) bias of f(xbar). It carries
    a Monte Carlo error of about (bootstrap error)/sqrt(numb_samples). The error is computed from the spread
    of func over the bootstrap samples. The output of func may be a scalar or a numpy object. For multidimensional 
    data, conf_axis says which axis holds the configurations to be resampled.

    Args:
        func (callable): Function that calculates the observable.
        data (array-like): Input data.
        numb_samples (int): Number of bootstrap samples.
        sample_size (int, optional): Size of each sample. Defaults to None, which uses the number of
          configurations.
        same_rand_for_obs (bool, optional): Every index before conf_axis labels an observable. If True, all
          observables are resampled with the same random configurations, i.e. they are treated as perfectly
          correlated. If False, each observable gets its own random configurations. Indices after conf_axis
          always share the random configuration. Defaults to False.
        conf_axis (int, optional): Axis to resample. Defaults to 1; it is set to 0 for one-dimensional data.
        return_sample (bool, optional): Also return the results from the individual samples? Defaults to False.
        seed (int, optional): Seed for the random generator. Defaults to None, which chooses a random seed.
        err_by_dist (bool, optional): Take the error from the 68% quantiles instead of the standard deviation?
          The error is then the distance from the median of the bootstrap distribution to its quantiles, i.e.
          a width of the distribution, applied around the bias-corrected central value. Defaults to True.
        args (tuple or dict, optional): Extra arguments for func. A dict is passed as **args. Defaults to ().
        nproc (int, optional): Number of threads. nproc=1 turns off parallelization. Defaults to DEFAULTTHREADS.

    Returns:
        samples (optionally), bias-corrected central value, bootstrap error
    """
    bts = nimbleBoot(func, data, numb_samples, sample_size, same_rand_for_obs, conf_axis, return_sample, seed,
                     err_by_dist, args, nproc)
    return bts.getResults()


class nimbleGaussianBoot:

    def __init__(self, func, data, data_std_dev, numb_samples, sample_size, same_rand_for_obs, return_sample, seed,
                 err_by_dist, useCovariance, Covariance, args, nproc, asym_err):

        checkType('int',numb_samples=numb_samples)
        checkType('int',sample_size=sample_size)
        checkType(bool,same_rand_for_obs=same_rand_for_obs)
        checkType(bool,return_sample=return_sample)
        checkType(bool,err_by_dist=err_by_dist)
        checkType(bool,useCovariance=useCovariance)
        checkType('int',nproc=nproc)
        checkType(bool,asym_err=asym_err)
        self._func=func
        self._data=np.array(data)
        self._data_std_dev=np.array(data_std_dev)
        self._numb_samples=numb_samples
        self._sample_size=sample_size
        self._same_rand_for_obs=same_rand_for_obs
        self._return_sample=return_sample
        self._seed=_autoSeed(seed)
        checkType('int',seed=self._seed)
        self._err_by_dist=err_by_dist
        self._useCovariance=useCovariance
        self._Covariance=Covariance
        self._args=args
        self._numb_observe = len(data)
        self._nproc = nproc

        # Factor F of each covariance, with F^T F = cov, computed as in rng.multivariate_normal
        if self._useCovariance:
            self._covFactors = []
            for k in range(self._numb_observe):
                if self._Covariance is None:
                    cov = np.diag(self._data_std_dev[k]**2)
                else:
                    cov = np.asarray(self._Covariance[k])
                _, s, vh = np.linalg.svd(cov)
                self._covFactors.append(np.sqrt(s)[:,None]*vh)

        self._sampleval=parallel_function_eval(self.getGaussianBootstrapEstimator,range(self._numb_samples),nproc=self._nproc,args=(self._seed,))

        # func of the means, i.e. of a sample with zero noise
        central_data = []
        for k in range(self._numb_observe):
            central_data.append(self._addNoise(k,np.zeros(self._sampleShape())))
        fbar = meanArgWrapper(self._func, np.array(central_data), self._args)
        self._mean = _biasCorrect(fbar, self._sampleval)
        if not self._err_by_dist:
            self._error = std_dev(self._sampleval)
        else:
            self._error = dev_by_dist(self._sampleval, return_both_q=asym_err)

    def __repr__(self) -> str:
        return "nimbleGaussianBoot"

    def _sampleShape(self):
        """
        Shape of one observable's sample.
        """
        if self._useCovariance:
            shape = (len(self._data[0]),)
            if self._sample_size != 1:
                shape = (self._sample_size,) + shape
        else:
            if self._sample_size == 1:
                shape = np.broadcast_shapes(np.shape(self._data[0]),np.shape(self._data_std_dev[0]))
            else:
                shape = (self._sample_size,)
        return shape

    def _addNoise(self,k,z):
        """
        Turn standard normals z into a sample of observable k.
        """
        if not self._useCovariance:
            return self._data[k] + self._data_std_dev[k]*z
        m = len(self._data[k])
        x = np.dot(z.reshape(-1,m), self._covFactors[k])
        x += self._data[k]
        return x.reshape(z.shape)

    def getGaussianBootstrapEstimator(self,i,my_seed):
        sample_data = []
        if self._same_rand_for_obs:
            # One set of random numbers, shared by all observables (perfectly correlated)
            z = TBRNG([my_seed,i]).standard_normal(self._sampleShape())
            for k in range(self._numb_observe):
                sample_data.append(self._addNoise(k,z))
        else:
            # Independent random numbers for each observable
            for k in range(self._numb_observe):
                z = TBRNG([my_seed,i,k]).standard_normal(self._sampleShape())
                sample_data.append(self._addNoise(k,z))

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
    """
    Gaussian (parametric) bootstrap. Like bootstr, but each sample is drawn from a normal distribution around
    the mean values in data, with width data_std_dev or covariance Covariance. The central value is the
    bias-corrected estimate 2 f(xbar) - <f*>, as in bootstr, with f(xbar) being func of the mean values.

    Args:
        func (callable): Function that calculates the observable.
        data (array-like): Mean value of each observable.
        data_std_dev (array-like): Standard deviation of each observable.
        numb_samples (int): Number of bootstrap samples.
        sample_size (int, optional): Number of draws per observable in each sample. If sample_size > 1, func
          has to average over the draws, and data_std_dev should be the standard deviation of a single
          measurement, not of the mean. Defaults to 1.
        same_rand_for_obs (bool, optional): Use the same random numbers for every observable, i.e. treat them
          as perfectly correlated? Defaults to False.
        return_sample (bool, optional): Also return the results from the individual samples? Defaults to False.
        seed (int, optional): Seed for the random generator. Defaults to None, which chooses a random seed.
        err_by_dist (bool, optional): Take the error from the 68% quantiles instead of the standard deviation?
          The error is then the distance from the median of the bootstrap distribution to its quantiles, i.e.
          a width of the distribution, applied around the bias-corrected central value. Defaults to True.
        useCovariance (bool, optional): Draw each observable from a multivariate normal? Defaults to False.
        Covariance (array-like, optional): Covariance matrix of each observable, used if useCovariance. Defaults
          to None, which uses diag(data_std_dev**2).
        args (tuple or dict, optional): Extra arguments for func. A dict is passed as **args. Defaults to ().
        nproc (int, optional): Number of threads. nproc=1 turns off parallelization. Defaults to DEFAULTTHREADS.
        asym_err (bool, optional): Return both distances from the median to the 68% quantiles, if err_by_dist,
          instead of the larger one. Defaults to False.

    Returns:
        samples (optionally), bias-corrected central value, bootstrap error
    """
    data = np.asarray(data)
    data_std_dev = np.asarray(data_std_dev)

    bts_gauss = nimbleGaussianBoot(func=func, data=data, data_std_dev=data_std_dev, numb_samples=numb_samples, 
                                   sample_size=sample_size, same_rand_for_obs=same_rand_for_obs,return_sample=return_sample, 
                                   seed=seed, err_by_dist=err_by_dist, useCovariance=useCovariance, Covariance=Covariance, 
                                   args=args, nproc=nproc, asym_err=asym_err)
    return bts_gauss.getResults()
