#
# utilities.py
#
# D. Clarke
#
# Some utilities that you might use in any program.
#

from subprocess import run, PIPE
import numpy as np
import time, re, datetime, os
import latqcdtools.base.logger as logger


# For byte conversions
_bytePrefix = { "Q"  : 1024**10,
                "R"  : 1024**9,
                "Y"  : 1024**8,  
                "Z"  : 1024**7, 
                "E"  : 1024**6, 
                "P"  : 1024**5, 
                "T"  : 1024**4, 
                "G"  : 1024**3, 
                "M"  : 1024**2, 
                "k"  : 1024, 
                "1"  : 1,
                1    : 1 }


def _getPrefix(byteString):
    if byteString=="B": 
        prefix=1
    else:
        prefix=byteString[0]
    return prefix


def _convert(text):
    if text.isdigit():
        return int(text)
    else:
        return text.lower()


def _alphanum_key(key):
    """ Splits the string `key` at any point where there is one or more consecutive digits. 
    The regular expression `([0-9]+)` is used to match one or more digits. The parentheses `()` 
    capture the matched digits as separate elements. For example, if `key` were 
    `'abc123def456ghi'`, the resulting list would be `['abc', '123', 'def', '456', 'ghi']`. """
    return [_convert(c) for c in re.split('([0-9]+)', key)]


# ---------------------------------------------------------------------------------- MAKE INTERNAL FUNCTIONS MORE SMOOTH


def isArrayLike(obj) -> bool:
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


def isHigherDimensional(obj) -> bool:
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
    """ Remove outermost brackets of array-like object, if possible. """ 
    if not isArrayLike(obj):
        return obj
    N = len(obj)
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


def toNumpy(*args):
    result = ()
    for obj in args:
        if isArrayLike(obj):
            obj = np.array(obj)
        result += (obj,)
    return result


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


def cleanOutput(*args,label=None) -> str:
    """ This method takes a bunch of args and formats them automatically for output. The idea is
    that you can use this method to ensure that columns are well lined up.

    Args:
        *args: The numbers you want to output, separated by commas. 
        label (str, optional): Put label to the left of your output. Defaults to None.

    Returns:
        str: formatted output string 
    """
    data = ()
    form = ''
    if label is not None:
        if not isinstance(label,str):
            logger.TBError('label must be a string')
        form += '%'+str(len(label))+'s'
        data += (label,)
    spacing = ''
    for col in args:
        if col is None:
            data += ('',)
            form += spacing+'%15s'
        elif isinstance(col,str):
            data += (col,)
            form += spacing+'%15s'
        elif isinstance(col,complex):
            data += (col.real,)
            data += (col.imag,)
            form += spacing+'%15.8e  %15.8e'
        elif isinstance(col,list):
            logger.TBError('Expected list of scalars rather than list of lists.')
        else:
            data += (col,)
            form += spacing+'%15.8e'
        spacing = '  '
    return form % data


def printClean(*args,label=None):
    """ Wrapper for cleanOutput that prints to screen.

    Args:
        *args: The numbers you want to output, separated by commas. 
        label (str, optional): Put label to the left of your output. Defaults to None.
    """
    logger.info(cleanOutput(*args,label))


def printDict(dic):
    """ Prints key, value pairs line by line. """
    if not isinstance(dic,dict):
        logger.TBError('Expected type', dict, 'but received', type(dic))
    for key in dic:
        if type(dic[key])==dict:
            logger.info(key)
            printDict(dic[key])
        else:
            printClean(key,dic[key])


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


def comesBefore(date1,date2,format="%Y/%m/%d %H:%M:%S") -> bool:
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


def naturalSort(l) -> list:
    """Sort list of strings so that, e.g. '10' comes after '9' rather than before it."""
    return sorted(l, key=_alphanum_key)


def find_nearest_idx(array, value) -> int:
    """ Find the index of the element of array nearest to value. """
    array = np.array(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def substringBetween(string,a,b) -> str:
    """ Find the substring of string between a and b. If a==b, it looks between the
    first and second occurences of a. 

    Args:
        string (str)
        a (str): starting delimiter 
        b (str): ending delimiter

    Returns:
        str: substring
    """
    start_index = string.index(a) + len(a)
    end_index   = string[start_index:].index(b) + start_index
    return string[start_index:end_index]


def deleteFile(target):
    """ Delete the file at target, if it exists. """
    if os.path.isfile(target):
        try:
            os.remove(target)
            logger.info("Deleted file",target)
            return
        except OSError:
            pass
    else:
        pass
    logger.warn('Unable to remove file',target)


def createFilePath(fullFileName):
    """ Create the directory path if it isn't there already. """
    if '/' in fullFileName:
        dir_path = os.path.dirname(fullFileName)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path,exist_ok=True)


def byteConvert(x,b1,b2):
    """ Convert between bytes given scientific prefixes.

    Args:
        x (float): Bytes in original units. 
        b1 (str): Original units. 
        b2 (str): Target units.

    Returns:
        float: Bytes in target units. 
    """
    p1 =_getPrefix(b1)
    p2 =_getPrefix(b2)
    num=_bytePrefix[p1]
    den=_bytePrefix[p2]
    return x*num/den


class timer:

    """ A class to facilitate doing rudimentary timings in the Toolbox. """

    def __init__(self):
        logger.info("Timer initialized.")
        self._tstart = time.time()
        self._tend   = self._tstart

    def __repr__(self) -> str:
        return "timer"

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