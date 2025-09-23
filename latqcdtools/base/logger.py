# 
# logger.py                                                               
# 
# H. Sandmeyer, D. Clarke 
# 
# Methods for logging and output.
# 

import sys, inspect, datetime, logging
from colorama import Fore


_PASS         = Fore.GREEN
_WARNING      = Fore.YELLOW 
_FAIL         = Fore.RED 
_ENDC         = '\033[0m'
RECORDLOG     = False
CURRENT_LEVEL = 4


class ToolboxException(Exception): pass

# These modules are not allowed to use the logger
_BANLIST = ['matplotlib']


_log_levels = {
    'ALL'      : 0,
    'DEBUG'    : 1,
    'DETAILS'  : 2,
    'PROGRESS' : 3,
    'INFO'     : 4,
    'WARN'     : 5,
    'NONE'     : 6
    }


def createLogFile(filename="Toolbox.log"):
    """ 
    Have output sent also to a log file filename. If this file already exists, it will get deleted. We use the
    logging module because it knows how to handle multiple processes writing to the same file. 
    """
    global RECORDLOG
    RECORDLOG = True
    # I recommend against changing the log level here. Other modules, such as Numba, also use this logger, and
    # hence you can get spammed if you set it to DEBUG. Similarly, keep the format the same, which otherwise will
    # print additional, unwanted strings at the beginning of each line.
    logging.basicConfig(filename=filename, encoding='utf-8', level=logging.INFO, format='%(message)s', filemode='w')
    for module in _BANLIST:
        logging.getLogger(module).setLevel(logging.CRITICAL)
    info('Created log file',filename)


def _log(outString):
    global RECORDLOG
    if RECORDLOG:
        logging.info(outString[:-1])


def _getCallerName(frame):
    """ 
    Gets the name of the function that calls the present function. 
    """
    # The way the frame works: Each nested function is labelled by the first index. 0 is the current function, i.e.
    # getCallerName. 1 is the function that called this function, and so on. Second index = 3 retrieves the name.
    currframe  = inspect.currentframe()
    callframe  = inspect.getouterframes(currframe, 2)
    methodName = str(callframe[frame][3])+': '
    # If the main script is calling it, we don't need a name.
    if methodName=='<module>: ':
        methodName='' 
    # If the method was called from within a class, get the class name too.
    call_locals = currframe.f_back.f_back.f_locals
    try:
        className = str(call_locals['self'])+'.'
    except KeyError:
        try:
            className = str(call_locals['locals_dict']['self'])+'.'
        except KeyError:
            className = ''
    return className+methodName


def _getTimeStamp():
    """ 
    Get HH:MM:SS 
    """
    return ' ['+datetime.datetime.now().strftime("%H:%M:%S")+']'


# ----------------------------------------------------------------------------------------------- DEPENDENT ON LOG LEVEL


def set_log_level(level):
    global CURRENT_LEVEL
    CURRENT_LEVEL = _log_levels[level]


def debug(*args,frame=2):
    global CURRENT_LEVEL
    if CURRENT_LEVEL <= 1:
        args   = [str(s) for s in args]
        output = _getTimeStamp()+' DEBUG: '+_getCallerName(frame)+(' '.join(args))
        print(output)
        _log(output+'\n')


def details(*args):
    global CURRENT_LEVEL
    if CURRENT_LEVEL <= 2:
        args   = [str(s) for s in args]
        output = _getTimeStamp()+' DETAILS: '+(' '.join(args))
        print(output)
        _log(output+'\n')


def progress(*args):
    global CURRENT_LEVEL
    if CURRENT_LEVEL <= 3:
        args   = [str(s) for s in args]
        output = _getTimeStamp()+' PROGRESS: '+(' '.join(args))
        print(output)
        _log(output+'\n')


def info(*args):
    global CURRENT_LEVEL
    if CURRENT_LEVEL <= 4:
        args   = [str(s) for s in args]
        output = _getTimeStamp()+' INFO: '+(' '.join(args))
        print(output)
        _log(output+'\n')


def warn(*args,frame=2):
    global CURRENT_LEVEL
    if CURRENT_LEVEL <= 5:
        args    = [str(s) for s in args]
        output  = _getTimeStamp()+' WARNING: '+_getCallerName(frame)+(' '.join(args))
        outputc = _getTimeStamp()+_WARNING+' WARNING: '+_getCallerName(frame)+(' '.join(args))+_ENDC
        print(outputc)
        _log(output+'\n')


# --------------------------------------------------------------------------------------------- INDEPENDENT OF LOG LEVEL


def TBFail(*args):
    args    = [str(s) for s in args]
    output  = _getTimeStamp()+' FAIL: '+(' '.join(args))
    outputc = _getTimeStamp()+_FAIL+' FAIL: '+(' '.join(args))+_ENDC
    print(outputc)
    _log(output + '\n')


def TBError(*args,frame=2):
    """
    Print error message and exit with -1.

    Args:
        frame (int, optional): Controls the name of the caller. Defaults to method that called TBError.
    """
    args    = [str(s) for s in args]
    output  = _getTimeStamp()+' ERROR: '+_getCallerName(frame)+(' '.join(args))
    outputc = _getTimeStamp()+_FAIL+' ERROR: '+_getCallerName(frame)+(' '.join(args))+_ENDC
    print(outputc)
    _log(output + '\n')
    sys.exit(-1)


def TBRaise(*args,frame=2,exception=None):
    """
    Add custom text to exceptions, which can be useful for debugging sometimes. 

    Args:
        frame (int, optional): Controls the name of the caller. Defaults to method that called TBRaise.
        exception (Exception, optional): If None, raise generic ToolboxException. Otherwise raise exception.
    """
    args    = [str(s) for s in args]
    output  = _getCallerName(frame)+(' '.join(args))
    outputc = _WARNING+_getCallerName(frame)+(' '.join(args))+_ENDC
    _log(output + '\n')
    if exception is None:
        raise ToolboxException(output) 
    else:
        print(outputc)
        raise exception 


def TBPass(*args):
    args    = [str(s) for s in args]
    output  = _getTimeStamp()+' SUCCESS: '+(' '.join(args))
    outputc = _getTimeStamp()+_PASS+' SUCCESS: '+(' '.join(args))+_ENDC
    print(outputc)
    _log(output + '\n')

