#
# check.py
#
# D. Clarke
#
# Warning and error control, along with methods to check internal code consistency. 
#

import math, warnings
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import envector,isArrayLike
warnings.filterwarnings("ignore", category=np.ComplexWarning)


class DivideByZeroError(Exception): pass
class UnderflowError(Exception): pass
class InvalidValueError(Exception): pass
class IllegalArgumentError(ValueError): pass


CATCHUNDERFLOW = True 


def err_handler(err, flag):
    """ This method lets us control in detail how different types of errors are treated. """
    global CATCHUNDERFLOW
    if flag == 1:
        raise DivideByZeroError(err)
    elif flag == 2:
        raise OverflowError
    elif flag == 4:
        if CATCHUNDERFLOW: 
            raise UnderflowError(err)
        else:
            pass
    elif flag == 8:
        raise InvalidValueError(err)
    elif flag == 9:
        raise DivideByZeroError(err)
    else:
        logger.TBError('Encountered unknown exception',err,'with flag',flag)


np.seterrcall(err_handler)
np.seterr(all='call')


def ignoreUnderflow():
    """ Turn off underflow crashes. """
    global CATCHUNDERFLOW
    CATCHUNDERFLOW = False
    logger.warn("Underflow behavior set to pass.")


def checkType(obj, expectedType):
    """ Check the type of an object.

    Args:
        obj (obj)
        expectedType (type): what type do you expect? Also accepts "array". 
    """
    if expectedType=="array":
        if not isArrayLike(obj):
            logger.TBError('Expected array-like object but received',type(obj),frame=3)
    else:
        if not isinstance(obj,expectedType):
            logger.TBError('Expected type',expectedType,'but received',type(obj),frame=3)


def checkDomain(obj, expectedDomain, variableName):
    """ Check that obj lies in expectedDomain.

    Args:
        obj (obj)
        expectedDomain (list): list of values obj is allowed to take
        variableName (str): name of obj 
    """    
    if not obj in expectedDomain:
        logger.TBError('Expected',variableName,'to be one of',expectedDomain)


def checkEqualLengths(*args):
    length = len(envector(args[0]))
    for array in args:
        if array is not None:
            if len(envector(array)) != length:
                logger.TBError('Array length mismatch detected.',frame=3)