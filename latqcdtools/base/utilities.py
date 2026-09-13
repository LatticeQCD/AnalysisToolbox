#
# utilities.py
#
# D. Clarke, K. Ebira
#
# Some utilities that you might use in any program.
#

from subprocess import run, PIPE
import numpy as np
import time, re, datetime
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
    try:
        # Try to convert text to a float for numerical sorting
        return float(text)
    except ValueError:
        # If it fails, return the text in lowercase
        return text.lower()


def _alphanum_key(key):
    """ 
    Splits the string `key` at any point where there's a change from digit to non-digit or vice versa.
    Will accurately handle integers and real numbers.
    """
    return [_convert(c) for c in re.split('([0-9]+(?:\\.[0-9]+)?)', key)]


# ---------------------------------------------------------------------------------- MAKE INTERNAL FUNCTIONS MORE SMOOTH


def isArrayLike(obj) -> bool:
    """ 
    Figure out whether obj is indexable.

    Args:
        obj (python object)

    Returns:
        bool: True if there is at least one index, false otherwise. 
    """
    if isinstance(obj, (dict, set)):
        return False
    try:
        obj[0]
        return True
    except (TypeError, IndexError, KeyError):
        return False


def isHigherDimensional(obj) -> bool:
    """ 
    Figure out whether obj has at least two indices.

    Args:
        obj (array-like)

    Returns:
        bool: True if there are at least two indices, false otherwise. 
    """
    if isinstance(obj, (dict, set)):
        return False
    try:
        obj[0][0]
        return True
    except (TypeError, IndexError, KeyError):
        return False


def isIntType(obj) -> bool:
    return isinstance(obj, (int, np.integer)) and not isinstance(obj, bool)


def isFloatType(obj) -> bool:
    return isinstance(obj, (float, np.floating))


def isComplexType(obj) -> bool:
    return isinstance(obj, (complex, np.complexfloating))


def isScalar(obj) -> bool:
    if (not isIntType(obj)) and (not isFloatType(obj)) and (not isComplexType(obj)):
        return False
    return True


def isReal(obj) -> bool:
    if not isScalar(obj): 
        return False
    if obj.imag != 0:
        return False
    return True


def unvector(obj):
    """ 
    Remove outermost brackets of array-like object with single element, if possible. This is needed
    because sometimes different numpy methods give inconsistent outputs, like turning a scalar
    into a zero-dimensional array, a 1-dimensional array, or just the scalar itself.

    Args:
        obj (python object)

    Returns:
        obj, obj[0], or obj.item() depending on obj
    
    """ 
    if isinstance(obj,np.ndarray):
        if obj.ndim==0:
            return obj.item()
    if not isArrayLike(obj):
        return obj
    N = len(obj)
    if N > 1:
        return obj
    else:
        return obj[0]


def envector(*args):
    """ 
    Change obj to a numpy array if it's a scalar. Sometimes required when, e.g., using np.vectorize. 
    """
    result = ()
    for obj in args:
        if not isArrayLike(obj):
            obj = np.array([obj])
        result += (obj,)
    return unvector(result)


def toNumpy(*args,**kwargs):
    result = ()
    for obj in args:
        if isArrayLike(obj):
            if len(obj)==0:
                logger.TBRaise('Received argument of length 0')
            obj = np.array(obj,**kwargs)
        result += (obj,)
    return result


def appendToDocstring(string=None,args=None,returns=None):
    def decorator(func):
        appendand=''
        if string is not None:
            appendand += string
        if args is not None:
            appendand += """\n    Args:""" + args
        if returns is not None:
            appendand += """\n    Returns:""" + returns
        if func.__doc__ is None:
            func.__doc__ = appendand 
        else:
            func.__doc__ = func.__doc__ + appendand 
        return func
    return decorator


# ------------------------------------------------------------------------------------------------- CONVENIENCE FOR USER


def getArgs(parser):
    """ 
    Get arguments from the ArgumentParser. Complain if you don't get exactly the correct arguments. 
    """
    args, invalid_args = parser.parse_known_args()
    if len(invalid_args)>0:
        logger.TBRaise("Received unrecognized arguments",invalid_args,".")
    return args


def printArg(message,param):
    """ 
    Some arguments are None by default, and you only want to print them if they are set. 
    """
    if param is not None:
        logger.info(message,param)


def cleanOutput(*args,label=None,sspace=20) -> str:
    """ 
    Format arguments automatically for lined-up column output.

    Parameters
    ----------
    *args : tuple
        Values to format (scalars or strings).
    label : str, optional
        Put label to the left of your output. Defaults to None.
    sspace : int or collections.abc.Iterable of int, optional
        Width for string columns. Can be an integer applied to all string
        columns, or an iterable of integers specifying per-string-column widths.
        Defaults to 20.

    Returns
    -------
    str
        Formatted output string.

    Raises
    ------
    ToolboxException
        If label is not a string, or if an array-like/nested object is passed in args,
        or if sspace is invalid or has insufficient entries.
    """
    if isIntType(sspace):
        sspace_list = None
        single_sspace = sspace
    elif not isinstance(sspace, str) and (isArrayLike(sspace) or hasattr(sspace, '__iter__')):
        sspace_list = list(sspace)
        for s in sspace_list:
            if not isIntType(s):
                logger.TBRaise('sspace elements must be integers.')
    else:
        logger.TBRaise('sspace must be an integer or iterable of integers.')

    data = ()
    form = ''
    if label is not None:
        if not isinstance(label,str):
            logger.TBRaise('label must be a string')
        form += '%'+str(len(label))+'s'
        data += (label,)
        spacing = '  '
    else:
        spacing = ''

    str_idx = 0
    for col in args:
        if col is None:
            data += ('',)
            form += spacing+'%15s'
        elif isinstance(col,str):
            if sspace_list is None:
                width = single_sspace
            else:
                if str_idx >= len(sspace_list):
                    logger.TBRaise('Not enough sspace entries for the number of string columns.')
                width = sspace_list[str_idx]
                str_idx += 1
            data += (col,)
            form += spacing+f'%{width}s'
        elif isComplexType(col):
            data += (col.real,)
            data += (col.imag,)
            form += spacing+'%15.8e  %15.8e'
        elif isArrayLike(col) or hasattr(col, '__iter__') or isinstance(col, (list, tuple, np.ndarray, set, dict, range)):
            logger.TBRaise('Expected scalars or strings rather than array-like objects.')
        elif isScalar(col):
            data += (col,)
            form += spacing+'%15.8e'
        else:
            logger.TBRaise('Expected scalars or strings rather than array-like objects.')
        spacing = '  '
    return form % data


def printClean(*args,label=None,sspace=20):
    """ 
    Wrapper for cleanOutput that prints formatted output to screen.

    Parameters
    ----------
    *args : tuple
        Values to format (scalars or strings).
    label : str, optional
        Put label to the left of your output. Defaults to None.
    sspace : int or collections.abc.Iterable of int, optional
        Width for string columns. Can be an integer applied to all string
        columns, or an iterable of integers specifying per-string-column widths.
        Defaults to 20.
    """
    logger.info(cleanOutput(*args,label=label,sspace=sspace).strip())


def printDict(dic,level=0):
    """ 
    Prints key, value pairs line by line. 
    """
    if not isinstance(dic,dict):
        logger.TBRaise('Expected type', dict, 'but received', type(dic))
    indent='  '*level
    for key in dic:
        if type(dic[key])==dict:
            logger.info(f'{indent}{key}:')
            printDict(dic[key],level+1)
        else:
            printClean(f'{indent}{key}:',dic[key])


def shell(*args):
    """ 
    Carry out the passed arguments args in the shell. Can be passed as a single
    string or as a list. Captures and returns output of shell command. E.g.
        shell('ls -lah')
    """
    args = [str(s) for s in args]
    command = ' '.join(args)
    process = run(command,shell=True,check=True,stdout=PIPE,universal_newlines=True)
    if process.returncode != 0:
        logger.TBRaise(f'{command} failed with return code {process.returncode}')
    return process.stdout


def shellVerbose(*args):
    """ 
    Same as shell, but instead of capturing output, print it to screen. 
    """
    args = [str(s) for s in args]
    command = ' '.join(args)
    process = run(command,shell=True,check=True,stdout=PIPE,universal_newlines=True)
    if process.returncode != 0:
        logger.TBRaise(f'{command} failed with return code {process.returncode}')
    print(process.stdout)


def comesBefore(date1,date2,format="%Y/%m/%d %H:%M:%S",beforeOrEqual=False) -> bool:
    """ 
    Check whether date1 comes before date2.

    Args:
        date1 (str)
        date2 (str)
        format (str): format for date strings. Defaults to "%Y/%m/%d %H:%M:%S"
        beforeOrEqual (bool): also return True if date1 == date2
        
    Returns:
        bool: date1 < date2 
    """
    date1_converted = datetime.datetime.strptime(date1, format)
    date2_converted = datetime.datetime.strptime(date2, format)
    if beforeOrEqual:
        return date1_converted <= date2_converted
    return date1_converted < date2_converted


def elapsedSeconds(date1,date2,format="%Y/%m/%d %H:%M:%S") -> float:
    """
    Compute elapsed time in seconds between date1 and date2. 

    Args:
        date1 (str)
        date2 (str)
        format (str): format for date strings. Defaults to "%Y/%m/%d %H:%M:%S"
        
    Returns:
        float: elapsed time in seconds 
    """
    date1_converted = datetime.datetime.strptime(date1, format)
    date2_converted = datetime.datetime.strptime(date2, format)
    elapsed = date1_converted-date2_converted
    return abs(elapsed.total_seconds())


def naturalSort(l) -> list:
    """
    Sort list of strings so that, e.g. '10' comes after '9' rather than before it.
    """
    return sorted(l, key=_alphanum_key)


def find_nearest_idx(array, value) -> int:
    """ 
    Find the index of the element of array nearest to value. 
    """
    array = np.array(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def substringBetween(string,a,b) -> str:
    """ 
    Find the substring of string between a and b. If a==b, it looks between the
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


def byteConvert(x,b1,b2):
    """ 
    Convert between bytes given scientific prefixes.

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

    """ 
    A class to facilitate doing rudimentary timings in the Toolbox. 
    """

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



# ------------------------------------------------------------------------------------------------- LEGACY WRAPPERS 

