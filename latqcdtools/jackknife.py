# 
# jackknife.py                                                               
# 
# H. Sandmeyer, D. Clarke 
# 
# This file contains a jackknife routine that can handle arbitrary return values of functions. 
# Please look first at jackknife_old to understand what is going on.
# 
import numpy as np
import math
from latqcdtools.statistics import std_mean, std_err, calc_cov
import latqcdtools.logger as logger
import concurrent.futures


''' Calculate pseudo values of elements of objects that are not tuple. '''
def pseudo(mean, mean_i, numb_blocks):
    return numb_blocks * mean - (numb_blocks - 1) * mean_i

def pseudo_val(mean, mean_i, numb_blocks):
    if type(mean) is tuple:
        retvalue = ()
        for i in range(len(mean)):
            retvalue += (pseudo(np.array(mean[i]), np.array(mean_i[i]), numb_blocks),)
        return retvalue
    else:
        return pseudo(np.array(mean), np.array(mean_i), numb_blocks)



class nimbleJack:

    """ Class allowing for parallelization of the jackknife function. """

    def __init__(self, func, data, nblocks, confAxis, return_sample, args, cov, parallelize, nproc):

        self._func=func
        self._data=np.array(data)
        self._nblocks=nblocks
        self._confAxis=confAxis
        self._return_sample=return_sample
        self._args=args
        self._cov=cov

        # If the measurements are accessed by the second index in data, we construct the jackknife 
        # manually to allow different size of sets of measurements. conf_axis==1 is the only case 
        # that allows different lengths of data arrays per observable.
        if self._confAxis == 1:
            try:
                self._lengths = [len(self._data[i]) for i in range(len(self._data))]
                self._blocksizes = [math.floor(length / self._nblocks) for length in self._lengths]
                for length in self._lengths:
                    if length < self._nblocks:
                        raise IndexError("More number of blocks than datapoints!")
            except TypeError:  # if we get an 1D array
                self._confAxis = 0
                self._length = self._data.shape[self._confAxis]
                self._blocksize = math.floor(self._length / self._nblocks)
        else:
            self._length = self._data.shape[self._confAxis]
            self._blocksize = math.floor(self._length / self._nblocks)

        if self._confAxis == 1:
            self._numb_observe = len(self._data)
            used_data = [self._data[i][:self._nblocks * self._blocksizes[i]] for i in range(self._numb_observe)]
            used_data = np.array(used_data)
        else:
            if self._length < self._nblocks:
                logger.TBError("More blocks than data. length, self._nblocks=", self._length, self._nblocks)
            rolled_data = np.rollaxis(self._data, self._confAxis)
            used_data = np.rollaxis(rolled_data[:self._blocksize*self._nblocks], 0, self._confAxis+1)

        if isinstance(args, dict):
            self._mean = func(used_data, **args)
        else:
            self._mean = func(used_data, *args)

        if parallelize:
            blockList=range(self._nblocks)
            with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
                blockval=executor.map(self.getJackknifeEstimator, blockList)
            self._blockval=list(blockval)
        else:
            blockval=[]
            for i in range(self._nblocks):
                blockval.append(self.getJackknifeEstimator(i))
            self._blockval=blockval

        self._mean = std_mean(self._blockval)
        if cov is False:
            self._error = std_err(self._blockval)
        else:
            self._cov = (1/len(self._blockval))*calc_cov(self._blockval)
            self._error = np.sqrt(self._cov.diagonal(axis1=-2,axis2=-1))


    ''' Gets the ith jackknife estimator from throwing away jackknife block i. '''
    def getJackknifeEstimator(self,i):
        block_data = []
        if self._confAxis == 1:
            for k in range(self._numb_observe):
                block_data.append(self._data[k][0:i*self._blocksizes[k]])
                block_data[k] = np.concatenate((block_data[k], self._data[k][(i + 1)
                    * self._blocksizes[k]:self._nblocks*self._blocksizes[k]]))
        else:
            # The Jackknife blocks are constructed by rolling the conf axis to the first index.
            # Then the jackknife blocks are built. Aferwards the we roll back
            rolled_data = np.rollaxis(self._data, self._confAxis)
            for j in range(i*self._blocksize):
                block_data.append(rolled_data[j])
            for j in range((i+1)*self._blocksize, self._nblocks*self._blocksize):
                block_data.append(rolled_data[j])
            block_data = np.rollaxis(np.array(block_data), 0, self._confAxis + 1)
        block_data = np.array(block_data)
        if isinstance(self._args, dict):
            mean_i = self._func(block_data, **self._args)
        else:
            mean_i = self._func(block_data, *self._args)
        mean_i = pseudo_val(self._mean, mean_i, self._nblocks)
        return mean_i


    def getResults(self):
      if self._cov is False:
          if self._return_sample:
              return self._blockval, self._mean, self._error
          else:
              return self._mean, self._error
      else:
          if self._return_sample:
              return self._blockval, self._mean, self._error, self._cov
          else:
              return self._mean, self._error, self._cov


def jackknife(func, data, numb_blocks=20, conf_axis=1, return_sample = False, args=(), cov=False, parallelize=True, nproc=32):
    """Jackknife routine for arbitray functions. This routine creates the jackknife like blocked
    subsets of data and passes them to the function in the same format as in the input data. So the 
    idea is to write a function that computes an observable from a given data set. This function 
    can be put into this jackkife routine and will get the jackknifed blocked data as input. Based 
    on the output of the function, the jackkife mean and error are computed. The function may return 
    multiple observables that are either scalars or numpy objects. You can pass a multidimensional 
    object as data, but the jackknife function has to know from which axis blocks should be removed. 
    This is controlled by conf_axis (default = 0 for one dimensional arrays and default = 1 for 
    higher order arrays). Look into the __init__ function of the nimbleJack class to see in detail
    how this is implemented.
        Parameters
        ----------
        func : callable
            The function that calculates the observable
    
        data : array_lik 
            Input data
    
        numb_blocks : integer
            Number of jackknife blocks
    
        conf_axis : integer, optional, default = 0 for dim(data) = 1
        and default = 1 for dim(data) >= 2
            Axis that should be resampled
   
        return_sample : boolean, optional, default = False                                           
        Along with the mean and the error also return the results from the individual samples 
        args : array_like or dict, default = ()
            optional arguements to be passed to func. If a dictionary the are passed as **args.

        parallelize : bool
            Should you make it parallel? You may need to turn it off for some reason.

        nproc : integer
            Number of threads to use if you choose to parallelize.
    """
    jk = nimbleJack(func, data, numb_blocks, conf_axis, return_sample, args, cov, parallelize, nproc)
    return jk.getResults()


def jackknife_old(func, data, return_sample=False, numb_blocks=20, args=()):
    """Older version of the jackknife that does not support objects, just scalars as a return 
    value of func. Furthermore it is not parallelized, so it is slower. However it is simpler 
    to read."""

    length = len(data[0])
    blocksize = math.floor(length / numb_blocks)

    ''' data is a 2D array. One dimension represents the observable being measured, while the 
        other dimension represents the number of measurements of each observable. '''
    numb_observe = len(data)

    ''' :X as an array argument means fetch everything up to, but not including, X. '''
    used_data = [data[i][:numb_blocks * blocksize] for i in range(numb_observe)]

    mean = func(used_data, *args)
    blockval = []
    for i in range(0, numb_blocks):
        block_data = []
        for k in range(0, numb_observe):
            block_data.append([])
            ''' These next two for loops construct a new block_data set, which 
                for each observable k includes all measurements except those falling
                into the ith jackknife block. Since this process continues up to
                only numb_blocks*blocksize, one sees that the last few measurements
                of the original data are sometimes omitted. '''
            for j in range(0, i * blocksize):
                block_data[k].append(data[k][j])
            for j in range((i + 1) * blocksize, numb_blocks * blocksize):
                block_data[k].append(data[k][j])
        ''' The function itself calculates a mean. We correct for the bias. '''
        mean_i = func(block_data, *args)
        mean_i = numb_blocks * mean - (numb_blocks - 1) * mean_i
        blockval.append(mean_i)
    mean = np.mean(blockval)
    error = 0.0
    for i in range(0, numb_blocks):
        error += (blockval[i] - mean) * (blockval[i] - mean)
    error = math.sqrt(error / (numb_blocks - 1))
    error /= math.sqrt(numb_blocks)
    if return_sample:
        return blockval, mean, error
    else:
        return mean, error


''' Demonstration how to use the jackknife function '''
def demonstrate_jackknife():
    def simple_mean(a):
        return np.mean(a)

    def div_old(a, b, c):
        return b * np.mean(a[0]) / (c * np.mean(a[1]))

    def f(a, b, c):
        return b * np.mean(a[0]) / (c * np.mean(a[1]))

    def div1(a, b, c):
        return ([[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]],
                [f(a, b, c), f(a, b, c)], f(a, b, c))

    def div2(a, b, c):
        return [[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]]

    def div3(a, b, c):
        return [f(a, b, c), f(a, b, c)]

    def div4(a, b, c):
        return f(a, b, c)

    def divnp1(a, b, c):
        return (np.array([[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]]),
                np.array([f(a, b, c), f(a, b, c)]), np.array(f(a, b, c)))

    def divnp2(a, b, c):
        return np.array([[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]])

    def divnp3(a, b, c):
        return np.array([f(a, b, c), f(a, b, c)])

    def divnp4(a, b, c):
        return np.array(f(a, b, c))

    def divnp5(a, b, c):
        return (np.array([[f(a, b, c), f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c), f(a, b, c)]]),
                np.array([f(a, b, c), f(a, b, c)]), np.array(f(a, b, c)))

    a, b, c = (np.random.normal(10, 2, size=1000), np.random.normal(10, 2, size=1000),
            np.random.normal(10, 2, size=1000))

    mn_simple, err_simple = jackknife(simple_mean, a, 10)
    print(mn_simple)
    print(err_simple)

    mn_old, err_old = jackknife_old(div_old, [a, b], 10, args=(2, 2))
    print(mn_old, err_old)

    mn1, err1 = jackknife(div1, [a, b], 10, args=(2, 2))
    print(mn1)
    print(err1)

    mn2, err2 = jackknife(div2, [a, b], 10, args=(2, 2))
    print(mn2)
    print(err2)

    mn3, err3 = jackknife(div3, [a, b], 10, args=(2, 2))
    print(mn3)
    print(err3)

    mn4, err4 = jackknife(div4, [a, b], 10, args=(2, 2))
    print(mn4)
    print(err4)

    mn1, err1 = jackknife(divnp1, [a, b], 10, args=(2, 2))
    print(mn1)
    print(err1)

    mn2, err2 = jackknife(divnp2, [a, b], 10, args=(2, 2))
    print(mn2)
    print(err2)

    mn3, err3 = jackknife(divnp3, [a, b], 10, args=(2, 2))
    print(mn3)
    print(err3)

    mn4, err4 = jackknife(divnp4, [a, b], 10, args=(2, 2))
    print(mn4)
    print(err4)

    mn4, err4 = jackknife(divnp5, [a, b], 10, args=(2, 2))
    print(mn4)
    print(err4)
