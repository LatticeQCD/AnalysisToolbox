# 
# logger.py                                                               
# 
# H. Sandmeyer, D. Clarke 
# 
# Methods for logging and output.
# 
import sys
from colorama import Fore


class bcolors:

    """ Colors for logging messages """

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


# ----------------------------------------------------------------------------------------------- DEPENDENT ON LOG LEVEL


def isLevel(level):
  return log_levels[level] >= current_level

def set_log_level(level):
  global current_level
  current_level = log_levels[level]

def log(level, *args, **kwargs):
  if log_levels[level] >= current_level:
    print(*args, **kwargs)

def debug(*args):
  if current_level <= 1:
    args = [str(s) for s in args]
    print('  DEBUG: '+(' '.join(args))) 

def details(*args):
  if current_level <= 2:
    args = [str(s) for s in args]
    print('  DETAILS: '+(' '.join(args))) 
        
def progress(*args):
  if current_level <= 3:
    args = [str(s) for s in args]
    print('  PROGRESS: '+(' '.join(args))) 

def info(*args):
  if current_level <= 4:
    args = [str(s) for s in args]
    print('  INFO: '+(' '.join(args))) 

def warn(*args):
  if current_level <= 5:
    args = [str(s) for s in args]
    print(bcolors.WARNING+'  WARNING: '+(' '.join(args))+bcolors.ENDC)


# --------------------------------------------------------------------------------------------- INDEPENDENT OF LOG LEVEL


def TBFail(*args):
  args = [str(s) for s in args]
  print(bcolors.FAIL+'  FAIL: '+(' '.join(args))+bcolors.ENDC) 

def TBError(*args):
  args = [str(s) for s in args]
  print(bcolors.FAIL+'  ERROR: '+(' '.join(args))+bcolors.ENDC) 
  sys.exit(-1)

def TBPass(*args):
  args = [str(s) for s in args]
  print(bcolors.PASS+'  SUCCESS: '+(' '.join(args))+bcolors.ENDC)

def TBRed(*args):
  args = [str(s) for s in args]
  print(bcolors.FAIL+(' '.join(args))+bcolors.ENDC) 

def TBGreen(*args):
  args = [str(s) for s in args]
  print(bcolors.PASS+(' '.join(args))+bcolors.ENDC)