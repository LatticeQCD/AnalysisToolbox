# 
# initialize.py                                                               
# 
# D. Clarke 
# 
# Some routines to set up the toolbox, especially for keeping a record of what you did. 
# 

import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import shell


INITIALIZED = False     # Global flag to check if initialization has already occurred.
DEFAULTSEED = 7271978   # Default seed for reproducibility (needed in testing). Do not Google this date.


def gitHash():
    """ Obtain the current git hash. This assumes the Toolbox has been correctly installed.

    Returns:
        str: git hash 
    """
    PYTHONPATH = shell('echo $PYTHONPATH')
    for entry in PYTHONPATH.split(':'):
        # The package was called AnalysisToolbox, then LatticeToolbox, then AnalysisToolbox again. 
        if ('LatticeToolbox' in entry) or ('AnalysisToolbox' in entry):
            toolboxLocation = entry.strip()
    hash=shell('git --git-dir="'+toolboxLocation+'/.git" rev-parse HEAD').strip()
    return hash 


def introduceYourself():
    """ Corporate branding. """
    logger.info()
    logger.info("  _          _   _   _         _____           _ _                ")
    logger.info(" | |    __ _| |_| |_(_) ___ __|_   _|__   ___ | | |__   _____  __ ")
    logger.info(" | |   / _` | __| __| |/ __/ _ \| |/ _ \ / _ \| | '_ \ / _ \ \/ / ")
    logger.info(" | |__| (_| | |_| |_| | (_|  __/| | (_) | (_) | | |_) | (_) >  <  ")
    logger.info(" |_____\__,_|\__|\__|_|\___\___||_|\___/ \___/|_|_.__/ \___/_/\_\ ")
    logger.info("                                                                  ")
    logger.info()


def initialize(logFile='Toolbox.log'):
    """ Some common tasks to do at the start of a run where you want to keep track of things. """
    global INITIALIZED
    INITIALIZED = True
    introduceYourself()
    logger.createLogFile(logFile)
    logger.info("Current git commit =",gitHash())
    logger.info()


def finalize():
    """ Some common tasks to do when you're done. """
    global INITIALIZED
    if not INITIALIZED:
        logger.warn('Called without having initialized first!')
    else:
        logger.info()
        logger.info("I'm finished!")
        logger.info()
