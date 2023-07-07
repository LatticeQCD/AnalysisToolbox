#
# check.py
#
# D. Clarke
#
# Warning and error control, along with methods to check internal code consistency. 
#

import warnings
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import envector, isArrayLike


# This warning is on by default, which worries that we may lose precision when using complex numbers. But
# I think in the instances I use complex numbers, they cannot really be avoided.
warnings.filterwarnings("ignore", category=np.ComplexWarning)


class DivideByZeroError(Exception): pass
class UnderflowError(Exception): pass
class InvalidValueError(Exception): pass
class IllegalArgumentError(ValueError): pass


CATCHUNDERFLOW    = True 
CATCHOVERFLOW     = True 
CATCHDIVIDEBYZERO = True
CATCHINVALIDVALUE = True


def err_handler(err, flag):
    """ This method lets us control in detail how different types of errors are treated. """
    global CATCHUNDERFLOW
    global CATCHOVERFLOW
    global CATCHDIVIDEBYZERO
    global CATCHINVALIDVALUE
    if flag == 1:
        if CATCHDIVIDEBYZERO:
            raise DivideByZeroError(err)
        else:
            pass
    elif flag == 2:
        if CATCHOVERFLOW:
            raise OverflowError
        else:
            pass
    elif flag == 4:
        if CATCHUNDERFLOW: 
            raise UnderflowError(err)
        else:
            pass
    elif flag == 8:
        if CATCHINVALIDVALUE:
            raise InvalidValueError(err)
        else:
            pass
    elif flag == 9:
        if CATCHDIVIDEBYZERO:
            raise DivideByZeroError(err)
        else:
            pass
    else:
        logger.TBError('Encountered unknown exception',err,'with flag',flag)


np.seterrcall(err_handler)
np.seterr(all='call')


def ignoreUnderflow():
    """ Turn off underflow crashes. """
    global CATCHUNDERFLOW
    CATCHUNDERFLOW = False
    logger.warn("Underflow behavior set to pass.")


def ignoreOverflow():
    """ Turn off overflow crashes. """
    global CATCHOVERFLOW
    CATCHOVERFLOW = False
    logger.warn("Overflow behavior set to pass.")


def ignoreDivideByZero():
    """ Turn off zero division crashes. """
    global CATCHDIVIDEBYZERO
    CATCHDIVIDEBYZERO = False
    logger.warn("Zero division behavior set to pass.")


def ignoreInvalidValue():
    """ Turn off invalid value crashes. """
    global CATCHINVALIDVALUE
    CATCHINVALIDVALUE = False
    logger.warn("Invalid value behavior set to pass.")


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
    """ Check that all array-like objects passed have the same length. """
    length = len(envector(args[0]))
    for array in args:
        if array is not None:
            if len(envector(array)) != length:
                logger.TBError('Array length mismatch detected.',frame=3)