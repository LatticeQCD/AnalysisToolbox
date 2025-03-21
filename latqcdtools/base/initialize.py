# 
# initialize.py                                                               
# 
# D. Clarke 
# 
# Some routines to set up the toolbox, especially for keeping a record of what you did. 
# 

import os, sys, warnings
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import shell, createFilePath


INITIALIZED = False     # Global flag to check if initialization has already occurred.
DEFAULTSEED = 7271978   # Default seed for reproducibility (needed in testing). Do not Google this date.


def gitHash() -> str:
    """ 
    Obtain the current git hash. This assumes the Toolbox has been correctly installed.

    Returns:
        str: git hash 
    """
    PYTHONPATH = shell('echo $PYTHONPATH')
    try:
        for entry in PYTHONPATH.split(':'):
            # The package was called AnalysisToolbox, then LatticeToolbox, then AnalysisToolbox again. 
            if ('AnalysisToolbox' in entry) or ('LatticeToolbox' in entry):
                toolboxLocation = entry.strip()
        hash=shell('git --git-dir="'+toolboxLocation+'/.git" rev-parse HEAD').strip()
    except:
        hash='GIT_NOT_FOUND'
    return hash 


def introduceYourself():
    """ 
    Corporate branding. ASCII generated from https://patorjk.com. 
    """
    logger.info()
    logger.info(r"     _                _           _    _____           _ _                ")
    logger.info(r"    / \   _ __   __ _| |_   _ ___(_)__|_   _|__   ___ | | |__   _____  __ ")
    logger.info(r"   / _ \ | '_ \ / _` | | | | / __| / __|| |/ _ \ / _ \| | '_ \ / _ \ \/ / ")
    logger.info(r"  / ___ \| | | | (_| | | |_| \__ \ \__ \| | (_) | (_) | | |_) | (_) >  <  ")
    logger.info(r" /_/   \_\_| |_|\__,_|_|\__, |___/_|___/|_|\___/ \___/|_|_.__/ \___/_/\_\ ")
    logger.info(r"                        |___/                                             ")
    logger.info()


def initialize(logFile=None):
    """ 
    Some common tasks to do at the start of a run where you want to keep track of things. 
    """
    global INITIALIZED
    warnings.simplefilter('ignore', UserWarning)
    introduceYourself()
    if logFile is None:
        logFile = 'log/' + os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.log'
    createFilePath(logFile)
    logger.createLogFile(logFile)
    logger.info("Current git commit =",gitHash())
    logger.info()
    INITIALIZED = True


def finalize():
    """ 
    Some common tasks to do when you're done. 
    """
    global INITIALIZED
    if not INITIALIZED:
        logger.warn('Called without having initialized first!')
    else:
        logger.info()
        logger.info("I'm finished!")
        logger.info()
