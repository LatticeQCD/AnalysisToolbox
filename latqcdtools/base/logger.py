# 
# logger.py                                                               
# 
# H. Sandmeyer, D. Clarke 
# 
# Methods for logging and output.
# 

import sys, inspect, datetime, logging
from colorama import Fore


PASS          = Fore.GREEN
WARNING       = Fore.YELLOW 
FAIL          = Fore.RED 
ENDC          = '\033[0m'
RECORDLOG     = False
CURRENT_LEVEL = 4


log_levels    = {
        'ALL' : 0,
        'DEBUG' : 1,
        'DETAILS' : 2,
        'PROGRESS' : 3,
        'INFO' : 4,
        'WARN' : 5,
        'NONE' : 6
        }


def createLogFile(filename="Toolbox.log"):
    """ Have output sent also to a log file filename. If this file already exists, it will get deleted. We use the
    logging module because it knows how to handle multiple processes writing to the same file. """
    global RECORDLOG
    RECORDLOG = True
    # I recommend against changing the log level here. Other modules, such as Numba, also use this logger, and
    # hence you can get spammed if you set it to DEBUG. Similarly, keep the format the same, which otherwise will
    # print additional, unwanted strings at the beginning of each line.
    logging.basicConfig(filename=filename, encoding='utf-8', level=logging.INFO, format='%(message)s', filemode='w')
    info('Created log file',filename)


def log(outString):
    global RECORDLOG
    if RECORDLOG:
        logging.info(outString[:-1])


def getCallerName(frame):
    """ Gets the name of the function that calls the present function. """
    currframe = inspect.currentframe()
    callframe = inspect.getouterframes(currframe, 2)
    # The way the frame works: Each nested function is labelled by the first index. 0 is the current function, i.e.
    # getCallerName. 1 is the function that called this function, and so on. Second index = 3 retrieves the name.
    return str(callframe[frame][3])


def getTimeStamp():
    """ Get HH:MM:SS """
    return ' ['+datetime.datetime.now().strftime("%H:%M:%S")+']'


# ----------------------------------------------------------------------------------------------- DEPENDENT ON LOG LEVEL


def isLevel(level):
    return log_levels[level] >= CURRENT_LEVEL


def set_log_level(level):
    global CURRENT_LEVEL
    CURRENT_LEVEL = log_levels[level]


def debug(*args,frame=2):
    if CURRENT_LEVEL <= 1:
        args       = [str(s) for s in args]
        callerName = getCallerName(frame)
        output     = getTimeStamp()+' DEBUG: '+callerName+'--'+(' '.join(args))
        print(output)
        log(output+'\n')


def details(*args):
    if CURRENT_LEVEL <= 2:
        args   = [str(s) for s in args]
        output = getTimeStamp()+' DETAILS: '+(' '.join(args))
        print(output)
        log(output+'\n')


def progress(*args):
    if CURRENT_LEVEL <= 3:
        args   = [str(s) for s in args]
        output = getTimeStamp()+' PROGRESS: '+(' '.join(args))
        print(output)
        log(output+'\n')


def info(*args):
    if CURRENT_LEVEL <= 4:
        args   = [str(s) for s in args]
        output = getTimeStamp()+' INFO: '+(' '.join(args))
        print(output)
        log(output+'\n')


def warn(*args,frame=2):
    if CURRENT_LEVEL <= 5:
        args       = [str(s) for s in args]
        callerName = getCallerName(frame)
        output     = getTimeStamp()+WARNING+' WARNING: '+callerName+'--'+(' '.join(args))+ENDC
        print(output)
        log(output+'\n')


# --------------------------------------------------------------------------------------------- INDEPENDENT OF LOG LEVEL


def TBFail(*args):
    args   = [str(s) for s in args]
    output = getTimeStamp()+FAIL+' FAIL: '+(' '.join(args))+ENDC
    print(output)
    log(output + '\n')


def TBError(*args,frame=2):
    args = [str(s) for s in args]
    callerName = getCallerName(frame)
    output = getTimeStamp()+FAIL+' ERROR: '+callerName+'--'+(' '.join(args))+ENDC
    print(output)
    log(output + '\n')
    sys.exit(-1)


def TBPass(*args):
    args = [str(s) for s in args]
    output = getTimeStamp()+PASS+' SUCCESS: '+(' '.join(args))+ENDC
    print(output)
    log(output + '\n')