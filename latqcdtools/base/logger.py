# 
# logger.py                                                               
# 
# H. Sandmeyer, D. Clarke 
# 
# Methods for logging and output.

# 

import sys, inspect, datetime, subprocess
from colorama import Fore


class bcolors:
    """ Colors for logging messages. """
    def __init__(self):
        pass
    PASS    = Fore.GREEN
    WARNING = Fore.YELLOW 
    FAIL    = Fore.RED 
    ENDC    = '\033[0m'


# Set the log level here.
current_level = 4
log_levels = {
        'ALL' : 0,
        'DEBUG' : 1,
        'DETAILS' : 2,
        'PROGRESS' : 3,
        'INFO' : 4,
        'WARN' : 5,
        'NONE' : 6
        }


recordLog = False
logFile   = None
def createLogFile(filename="Toolbox.log"):
    """ Have output sent also to a log file. """
    global recordLog
    global logFile
    logFile = open(filename,'w')
    recordLog = True


def closeLogFile():
    """ Close the log file. """
    global recordLog
    global logFile
    if not recordLog:
        warn('Attempted to close log file with makeLog=False.')
    logFile.close()
    recordLog = False


def log(outString):
    global recordLog
    global logFile
    if recordLog:
        logFile.write(outString)


def getCallerName(frame):
    """ Gets the name of the function that calls the present function. """
    currframe = inspect.currentframe()
    callframe = inspect.getouterframes(currframe, 2)
    # The way the frame works: Each nested function is labelled by the first index. 0 is the current function, i.e.
    # getCallerName. 1 is the function that called this function, and so on. Second index = 3 retrieves the name.
    callerName = str(callframe[frame][3])
    return callerName


def getTimeStamp():
    """ Get HH:MM:SS """
    return ' ['+datetime.datetime.now().strftime("%H:%M:%S")+']'


# ----------------------------------------------------------------------------------------------- DEPENDENT ON LOG LEVEL


def isLevel(level):
    return log_levels[level] >= current_level


def set_log_level(level):
    global current_level
    current_level = log_levels[level]


def debug(*args,frame=2):
    if current_level <= 1:
        args       = [str(s) for s in args]
        callerName = getCallerName(frame)
        output     = getTimeStamp()+' DEBUG: '+callerName+'--'+(' '.join(args))
        print(output)
        log(output+'\n')


def details(*args):
    if current_level <= 2:
        args   = [str(s) for s in args]
        output = getTimeStamp()+' DETAILS: '+(' '.join(args))
        print(output)
        log(output+'\n')


def progress(*args):
    if current_level <= 3:
        args   = [str(s) for s in args]
        output = getTimeStamp()+' PROGRESS: '+(' '.join(args))
        print(output)
        log(output+'\n')


def info(*args):
    if current_level <= 4:
        args   = [str(s) for s in args]
        output = getTimeStamp()+' INFO: '+(' '.join(args))
        print(output)
        log(output+'\n')


def warn(*args,frame=2):
    if current_level <= 5:
        args       = [str(s) for s in args]
        callerName = getCallerName(frame)
        output     = getTimeStamp()+bcolors.WARNING+' WARNING: '+callerName+'--'+(' '.join(args))+bcolors.ENDC
        print(output)
        log(output+'\n')


# --------------------------------------------------------------------------------------------- INDEPENDENT OF LOG LEVEL


def TBFail(*args):
    args   = [str(s) for s in args]
    output = getTimeStamp()+bcolors.FAIL+' FAIL: '+(' '.join(args))+bcolors.ENDC
    print(output)
    log(output + '\n')


def TBError(*args,frame=2):
    args = [str(s) for s in args]
    callerName = getCallerName(frame)
    output = getTimeStamp()+bcolors.FAIL+' ERROR: '+callerName+'--'+(' '.join(args))+bcolors.ENDC
    print(output)
    log(output + '\n')
    sys.exit(-1)


def TBPass(*args):
    args = [str(s) for s in args]
    output = getTimeStamp()+bcolors.PASS+' SUCCESS: '+(' '.join(args))+bcolors.ENDC
    print(output)
    log(output + '\n')


def gitHash():
    """ Obtain the current git hash.

    Returns:
        str: git hash 
    """
    process = subprocess.Popen(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE)
    output, _ = process.communicate()
    return output.decode('UTF-8').strip()


def introduceYourself():
    """ Print name along with git hash. """
    info()
    info("                       _           _  _______          _ _                ") 
    info("     /\               | |         (_)|__   __|        | | |               ") 
    info("    /  \   _ __   __ _| |_   _ ___ _ ___| | ___   ___ | | |__   _____  __ ") 
    info("   / /\ \ | '_ \ / _` | | | | / __| / __| |/ _ \ / _ \| | '_ \ / _ \ \/ / ") 
    info("  / ____ \| | | | (_| | | |_| \__ \ \__ \ | (_) | (_) | | |_) | (_) >  <  ") 
    info(" /_/    \_\_| |_|\__,_|_|\__, |___/_|___/_|\___/ \___/|_|_.__/ \___/_/\_\ ") 
    info("                          __/ |                                           ") 
    info("                         |___/                                            ") 
    info()
    info("Current git commit =",gitHash())
    info()
