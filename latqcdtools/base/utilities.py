#
# utilities.py
#
# D. Clarke
#
# Some utilities that you might use in any program.
#

from subprocess import run, PIPE
import numpy as np
import time, re, datetime, os, shutil
import latqcdtools.base.logger as logger
import glob


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
    """ 
    Splits the string `key` at any point where there is one or more consecutive digits. 
    The regular expression `([0-9]+)` is used to match one or more digits. The parentheses `()` 
    capture the matched digits as separate elements. For example, if `key` were 
    `'abc123def456ghi'`, the resulting list would be `['abc', '123', 'def', '456', 'ghi']`. 
    """
    return [_convert(c) for c in re.split('([0-9]+)', key)]


# ---------------------------------------------------------------------------------- MAKE INTERNAL FUNCTIONS MORE SMOOTH


def isArrayLike(obj) -> bool:
    """ 
    Figure out whether obj is indexable.

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
    """ 
    Figure out whether obj has at least two indices.

    Args:
        obj (array-like)

    Returns:
        bool: True if there are at least two indices, false otherwise. 
    """
    try:
        obj[0][0]
        return True
    except (TypeError, IndexError):
        return False


def isIntType(obj) -> bool:
    return isinstance(obj,(int, np.int8, np.int16, np.int32, np.int64))


def isFloatType(obj) -> bool:
    return isinstance(obj,(float,np.float16,np.float32,np.float64,np.float128))


def isComplexType(obj) -> bool:
    return isinstance(obj,(complex,np.complex64,np.complex128))


def isScalar(obj) -> bool:
    if (not isIntType(obj)) and (not isFloatType(obj)) and (not isComplexType(obj)):
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


def cleanOutput(*args,label=None) -> str:
    """ 
    This method takes a bunch of args and formats them automatically for output. The idea is
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
            logger.TBRaise('label must be a string')
        form += '%'+str(len(label))+'s'
        data += (label,)
    spacing = ''
    for col in args:
        if col is None:
            data += ('',)
            form += spacing+'%15s'
        elif isinstance(col,str):
            data += (col,)
            form += spacing+'%20s'
        elif isinstance(col,complex):
            data += (col.real,)
            data += (col.imag,)
            form += spacing+'%15.8e  %15.8e'
        elif isinstance(col,list):
            logger.TBRaise('Expected list of scalars rather than list of lists.')
        else:
            data += (col,)
            form += spacing+'%15.8e'
        spacing = '  '
    return form % data


def printClean(*args,label=None):
    """ 
    Wrapper for cleanOutput that prints to screen.

    Args:
        *args: The numbers you want to output, separated by commas. 
        label (str, optional): Put label to the left of your output. Defaults to None.
    """
    logger.info(cleanOutput(*args,label).strip())


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
    process = run(' '.join(args),shell=True,check=True,stdout=PIPE,universal_newlines=True)
    return process.stdout


def shellVerbose(*args):
    """ 
    Same as shell, but instead of capturing output, print it to screen. 
    """
    args = [str(s) for s in args]
    process = run(' '.join(args),shell=True,check=True,stdout=PIPE,universal_newlines=True)
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


def deleteLine(target,line_number):
    """
    Delete line line_number from file target.

    Args:
        target (str)
        line_number (int)
    """
    if os.path.isfile(target):
        with open(target, 'r') as file:
            lines = file.readlines()
        if 0 < line_number <= len(lines):
            lines.pop(line_number - 1)  # line_number is 1-based
            with open(target, 'w') as file:
                file.writelines(lines)
        else:
            logger.TBRaise("Line number is out of range.")
        return
    logger.warn(f"{target} does not exist.")


def deleteFile(target):
    """ 
    Delete the file at target, if it exists. 
    """
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


def deleteFolder(target):
    """ 
    Delete the folder at target, if it exists. 
    """
    if os.path.isdir(target):
        try:
            shutil.rmtree(target)
            logger.info("Deleted folder",target,"and its subdirectories.")
            return
        except OSError:
            pass
    else:
        pass
    logger.warn('Unable to remove folder',target)


def createFilePath(fullFileName):
    """ 
    Create the directory path if it isn't there already. 
    """
    if '/' in fullFileName:
        dir_path = os.path.dirname(fullFileName)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path,exist_ok=True)


def ls(filePath) -> list:
    return naturalSort(list(glob.iglob(filePath)))


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