# 
# speedify.py                                                               
# 
# L. Altenkort, D. Clarke 
# 
# Some methods and classes to easily make python code faster. 
#


import numpy as np
import concurrent.futures
from latqcdtools.base.utilities import envector, getMaxThreads
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger
from numba import njit
from numba.typed import List


COMPILENUMBA = False


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
        logger.details('Compiling',func.__name__+'.')
        return njit(func)
    else:
        return func


def ennumba(input):
    """ Convert a 1-d array-like input into a numba List object. Needed whenever you want to pass
    a 1-d array-like object to numba. I recommend to put this outside your numba function, because
    if you have to keep doing this conversion each time the function is called, you will pay for it.

    Args:
        input (array-like)

    Returns:
        numba.typed.List: List version of input array 
    """    
    output = List()
    [output.append(x) for x in input]
    return output    


class ComputationClass:

    """ A class to parallelize functions. To be used with parallel_function_eval for convenience. """

    def __init__(self, function, input_array, nproc, *add_param):
        checkType(nproc,int)
        checkType(input_array,"array")
        self._nproc = nproc  # number of processes
        self._input_array = input_array
        self._function = function
        if nproc > getMaxThreads():
            logger.warn('We recommend using fewer processes than',getMaxThreads()) 

        # additional arguments for actual_computation
        self._add_param = add_param

        # compute the result when class is initialized
        self._result = self.parallelization_wrapper()

    def parallelization_wrapper(self):
        results = []
        if self._nproc==1:
            for i in self._input_array:
                results.append(self.pass_argument_wrapper(i)) 
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
                for result in executor.map(self.pass_argument_wrapper, self._input_array):
                    results.append(list(envector(result)))
            results = list(map(list, zip(*results)))  # "transpose" the list to allow for multiple return values like a normal function.
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
        logger.info('Using for-loop instead of concurrent.futures.')
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