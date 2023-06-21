#
# utilities.py
#
# D. Clarke
#
# Some utilities that you might use in any program.
#

from subprocess import run, PIPE
import numpy as np
import time, re, datetime
import latqcdtools.base.logger as logger


# ---------------------------------------------------------------------------------- MAKE INTERNAL FUNCTIONS MORE SMOOTH


def convertToNumpy(*args):
    """ Change lists or tuples to numpy arrays.

    Returns:
        tuple: tuple of args, each arg converted to a numpy array. 
    """
    args=list(args)
    for i in range(len(args)):
        if isinstance(args[i], (list, tuple)):
            args[i] = np.array(args[i])
    return tuple(args)


def isArrayLike(obj):
    """ Figure out whether obj is indexable.

    Args:
        obj (python object)

    Returns:
        bool: True if there is at least one index, false otherwise. 
    """
    try:
        obj[0]
        return True
    except (TypeError,IndexError):
        return False


def isHigherDimensional(obj):
    """ Figure out whether obj has at least two indices.

    Args:
        data (array-like)

    Returns:
        bool: True if there are at least two indices, false otherwise. 
    """
    try:
        obj[0][0]
        return True
    except (TypeError, IndexError):
        return False


def unvector(obj):
    """ Change obj to a scalar if it's an array-like object of length 1. Otherwise don't do anything. """
    try:
        N = len(obj)
    except TypeError:
        return obj
    if N > 1:
        return obj
    else:
        return obj[0]


def envector(*args):
    """ Change obj to a numpy array if it's a scalar. Sometimes required when, e.g., using np.vectorize. """
    result = ()
    for obj in args:
        if not isArrayLike(obj):
            obj = np.array([obj])
        result += (obj,)
    return unvector(result)


# ------------------------------------------------------------------------------------------------- CONVENIENCE FOR USER


def getArgs(parser):
    """ Get arguments from the ArgumentParser. Complain if you don't get exactly the correct arguments. """
    args, invalid_args = parser.parse_known_args()
    if len(invalid_args)>0:
        logger.TBError("Received unrecognized arguments",invalid_args,".")
    return args


def printArg(message,param):
    """ Some arguments are None by default, and you only want to print them if they are set. """
    if param is not None:
        logger.info(message,param)


def printDict(dic):
    """ Prints key, value pairs line by line. """
    if not isinstance(dic,dict):
        logger.TBError('Expected type', dict, 'but received', type(dic))
    for key in dic:
        logger.info(key,dic[key])


def printClean(*args,label=None):
    data = ()
    form = ''
    if label is not None:
        if not isinstance(label,str):
            logger.TBError('label must be a string')
        form += '%'+str(len(label))+'s: '
        data += (label,)
    for col in args:
        if col is None:
            data += ('',)
            form += '%16s'
        elif isinstance(col,complex):
            data += (col.real,)
            data += (col.imag,)
            form += '(%.8e  %.8e)  '
        else:
            data += (col,)
            form += '%.8e  '
    try:
        logger.info(form % data)
    except TypeError:
        logger.TBError('args',args,'is not list-type of floats.')


def shell(*args):
    """ Carry out the passed arguments args in the shell. Can be passed as a single
        string or as a list. Captures and returns output of shell command. E.g.
            shell('ls -lah')
    """
    args = [str(s) for s in args]
    process = run(' '.join(args),shell=True,check=True,stdout=PIPE,universal_newlines=True)
    return process.stdout


def shellVerbose(*args):
    """ Same as shell, but instead of capturing output, print it to screen. """
    args = [str(s) for s in args]
    process = run(' '.join(args),shell=True,check=True,stdout=PIPE,universal_newlines=True)
    print(process.stdout)


def getMaxThreads():
    """ Figure out how many threads are on this system. """
    return int( shell("lscpu | awk '/CPU\\(s\\)/ {print $2; exit}'") )


def comesBefore(date1,date2,format="%Y/%m/%d %H:%M:%S"):
    """ Check whether date1 comes before date2.

    Args:
        date1 (str)
        date2 (str)
        format (str): format for date strings. Defaults to "%Y/%m/%d %H:%M:%S"

    Returns:
        bool: date1 < date2 
    """
    date1_converted = datetime.datetime.strptime(date1, format)
    date2_converted = datetime.datetime.strptime(date2, format)
    return date1_converted < date2_converted


def naturalSort(l):
    """ Sort list of strings so that, e.g. '10' comes after '9' rather than before it. """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def find_nearest_idx(array, value):
    """ Find the index of the element of array nearest to value. """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class timer:

    """ A class to facilitate doing rudimentary timings in the Toolbox. """

    def __init__(self):
        logger.info("Timer initialized.")
        self._tstart = time.time()
        self._tend   = self._tstart


    def printTiming(self, message=None):
        self._tstart = self._tend
        self._tend   = time.time()
        timing = self._tend - self._tstart
        if message is None:
            logger.info("Time to finish: %12.8f [s]." % timing)
        else:
            logger.info("Time to finish "+message+": %12.8f [s]." % timing)


    def resetTimer(self):
        self._tstart = time.time()
        self._tend   = self._tstart