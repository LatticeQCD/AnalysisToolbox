# 
# speedify.py                                                               
# 
# L. Altenkort, D. Clarke 
# 
# Some methods and classes to easily make python code faster. 
#

import os
import numpy as np
import concurrent.futures
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger
from numba import njit
from numba.typed import List


# Resolve parallelizer dependencies
DEFAULTPARALLELIZER = 'pathos.pools'
try:
    import pathos.pools
except ModuleNotFoundError:
    DEFAULTPARALLELIZER = 'concurrent.futures'


COMPILENUMBA        = False
MAXTHREADS          = os.cpu_count() 
DEFAULTTHREADS      = MAXTHREADS - 2


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


class ComputationClass:

    def __init__(self, function, input_array, args, nproc, parallelizer):
        """ This class contains everything needed to parallelize a function. It allows for you to
        use your favorite parallelization library and pass arguments to the function.

        Args:
            function (func): to-be-parallelized function
            input_array (array-like): run function over this array 
            nproc (int): number of processors 
            parallelizer (str): Which library should I use to parallelize?
            *add_param: Pass any additional parameters as you would to the function
        """
        checkType(nproc,int)
        checkType(input_array,"array")
        checkType(parallelizer,str)
        self._input_array  = input_array
        self._function     = function
        self._parallelizer = parallelizer
        self._args         = args
        self._nproc        = nproc
        logger.debug('Initializing')
        logger.debug('input_array =',self._input_array)
        logger.debug('args =',self._args)
        if nproc > MAXTHREADS:
            logger.warn('We recommend using fewer processes than',MAXTHREADS) 
        self._result = self.parallelization_wrapper() # compute the result when class is initialized

    def __repr__(self) -> str:
        return "ComputationClass"

    def parallelization_wrapper(self):
        if self._nproc==1:
            results = []
            for i in self._input_array:
                results.append(self.pass_argument_wrapper(i)) 
        else:
            if self._parallelizer=='concurrent.futures':
                with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
                    results=executor.map(self.pass_argument_wrapper, self._input_array)
            elif self._parallelizer=='pathos.pools':
                pool = pathos.pools.ProcessPool(processes=self._nproc)
                results = pool.map(self.pass_argument_wrapper, self._input_array)
                pool.close()
                pool.join() 
                pool.clear()
            else:
                logger.TBError('Unknown parallelizer',self._parallelizer)
            results = list(results)
        return results

    def pass_argument_wrapper(self, single_input):
        return self._function(single_input, *self._args)

    def getResult(self):
        return self._result


def parallel_function_eval(function, input_array, args=(), nproc=DEFAULTTHREADS, parallelizer=DEFAULTPARALLELIZER):
    """ Parallelize a function over an input_array. Effectively this can replace a loop over an array and should
    lead to a performance boost.

    Args:
        function (func): to-be-parallelized function 
        input_array (array-like): array over which it should run 
        nproc (int): number of processes 

    Returns:
        array-like: func(input_array)
    """
    computer = ComputationClass(function, input_array, args=args, nproc=nproc, parallelizer=parallelizer)
    if nproc==1:
        logger.details('Using for-loop instead of',parallelizer)
    return computer.getResult()


def parallel_reduce(function, input_array, args=(), nproc=DEFAULTTHREADS, parallelizer=DEFAULTPARALLELIZER) -> float:
    """ Parallelize a function over an input_array, then sum over the input_array elements. 

    Args:
        function (func): to-be-parallelized function 
        input_array (array-like): array over which it should run 
        nproc (int): number of processes 

    Returns:
        float
    """
    container = parallel_function_eval(function, input_array, nproc=nproc, args=args, parallelizer=parallelizer)
    return np.sum(container)