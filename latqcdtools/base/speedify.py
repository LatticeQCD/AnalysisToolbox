# 
# speedify.py                                                               
# 
# L. Altenkort, D. Clarke 
# 
# Some methods and classes to easily make python code faster. 
#


import numpy as np
import concurrent.futures
import pathos.pools
from latqcdtools.base.utilities import shell 
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger
from numba import njit
from numba.typed import List


def getMaxThreads():
    """ Figure out how many threads are on this system. """
    return int( shell("lscpu | awk '/CPU\\(s\\)/ {print $2; exit}'") )


COMPILENUMBA = False
MAXTHREADS   = getMaxThreads()


def numbaON():
    """ Use numba wherever possible. By default it is turned off, since compilation takes some time,
    and hence you will only see a performance boost for particularly long-running functions. Must be
    called at the beginning of your code. """
    global COMPILENUMBA
    COMPILENUMBA = True
    logger.info('Using numba to speed things up.')


def numbaOFF():
    """ Turn off numba compilation for small functions. """ 
    global COMPILENUMBA
    COMPILENUMBA = False 
    logger.info('No longer using numba.')


def compile(func):
    global COMPILENUMBA
    if COMPILENUMBA:
        logger.info('Compiling',func.__name__+'.')
        return njit(func)
    else:
        return func


def numbaList(inList):
    """ Turn a list into List that numba can parse. """ 
    global COMPILENUMBA
    if COMPILENUMBA:
        nList = List()
        [nList.append(x) for x in inList]
        return nList 
    else:
        return inList


def setNproc(parallelize,nproc):
    """ If nproc = 1, then parallel_function_eval will do a for-loop behind the scenes. """
    if parallelize:
        return nproc
    else:
        return 1


class ComputationClass:

    def __init__(self, function, input_array, nproc, parallelizer='concurrent.futures', *add_param):
        checkType(nproc,int)
        checkType(input_array,"array")
        checkType(parallelizer,str)
        self._input_array  = input_array
        self._function     = function
        self._nproc        = nproc
        self._parallelizer = parallelizer
        if nproc > MAXTHREADS:
            logger.warn('We recommend using fewer processes than',MAXTHREADS) 

        self._add_param = add_param
        # compute the result when class is initialized
        self._result = self.parallelization_wrapper()

    def parallelization_wrapper(self):
        if self._nproc==1:
            results = []
            for i in self._input_array:
                results.append(self.pass_argument_wrapper(i)) 
        else:
#            with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
#                results=executor.map(self.pass_argument_wrapper, self._input_array)
#            results = list(results)
            with ProcessPool(processes=self._nproc) as pool:
                results = pool.map(self.pass_argument_wrapper, self._input_array)
            results = list(results)
        return results

    def pass_argument_wrapper(self, single_input):
        return self._function(single_input, *self._add_param)

    def getResult(self):
        return self._result


def parallel_function_eval(function, input_array, nproc, *add_param):
    """ Parallelize a function over an input_array. Effectively this can replace a loop over an array and should
    lead to a performance boost.

    Args:
        function (func): to-be-parallelized function 
        input_array (array-like): array over which it should run 
        nproc (int): number of processes 

    Returns:
        array-like: func(input_array)
    """
    computer = ComputationClass(function, input_array, nproc, *add_param)
    if nproc==1:
        logger.details('Using for-loop instead of concurrent.futures.')
    return computer.getResult()


def parallel_reduce(function, input_array, nproc, *add_param):
    """ Parallelize a function over an input_array, then sum over the input_array elements. 

    Args:
        function (func): to-be-parallelized function 
        input_array (array-like): array over which it should run 
        nproc (int): number of processes 

    Returns:
        float-like
    """
    container=parallel_function_eval(function, input_array, nproc, *add_param)
    return np.sum(container)