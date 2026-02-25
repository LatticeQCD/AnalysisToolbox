# 
# fileSystem                                                               
# 
# D. Clarke
# 
# Some bash-like operations. Sometimes it is better to just write a
# Python script to do these things.
#
import os, shutil, glob, datetime
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.base.utilities import naturalSort


def deleteLine(target,line_number):
    """
    Delete line line_number from file target, indexed from 1.

    Args:
        target (str)
        line_number (int)
    """
    checkType(str,target=target)
    checkType('int',line_number=line_number)
    if os.path.isfile(target):
        with open(target, 'r') as file:
            lines = file.readlines()
        if 0 < line_number <= len(lines):
            lines.pop(line_number - 1)  # line_number is indexed from 1
            with open(target, 'w') as file:
                file.writelines(lines)
        else:
            logger.TBRaise("Line number is out of range.")
        return
    logger.warn(f"{target} is not a regular file.")


def rm(target):
    """
    Delete target regular file or folder. Equivalent to rm -rf in bash.
    
    Args:
        target (str)
    """
    checkType(str,target=target)
    if os.path.isfile(target):
        try:
            os.remove(target)
            logger.info(f"Deleted regular file {target}")
            return
        except OSError:
            pass
    elif os.path.isdir(target):
        try:
            shutil.rmtree(target)
            logger.info(f"Deleted folder {target} and its subdirectories.")
            return
        except OSError:
            pass
    else:
        pass
    logger.warn(f"Unable to remove {target}")


def ls(target) -> list:
    """
    Get list of files in file path. Similar to ls in bash.

    Args:
        target (str)

    Returns:
        list: List of files in target 
    """
    checkType(str,target=target)
    return naturalSort(list(glob.iglob(target)))


def createFilePath(target):
    """ 
    Create the directory path for a file if it isn't there already. 

    Args:
        target (str)
    """
    checkType(str,target=target)
    if '/' in target:
        dir_path = os.path.dirname(target)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path,exist_ok=True)


def getFileTimeStamp(target) -> str:
    """
    Get the time stamp (when it was last modified) of a regular file.

    Args:
        target (str)

    Returns:
        str: time stamp in format 2025-09-23 14:56:27
    """
    checkType(str,target=target)
    if os.path.isfile(target):
        modification_time = os.path.getmtime(target)
        mod_time_readable = datetime.datetime.fromtimestamp(modification_time)
        return str(mod_time_readable)
    else:
        logger.warn(f"{target} is not regular file.")
    