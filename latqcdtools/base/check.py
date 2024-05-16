#
# check.py
#
# D. Clarke
#
# Warning and error control, along with methods to check internal code consistency. 
#

import warnings, inspect
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import envector, isArrayLike


# This warning is on by default, which worries that we may lose precision when using complex numbers. But
# I think in the instances I use complex numbers, they cannot really be avoided.
warnings.filterwarnings("ignore", category=np.ComplexWarning)

# I want solvers to try other strategies when they hit RuntimeWarnings.
warnings.filterwarnings("error", category=RuntimeWarning)


class DivideByZeroError(Exception): pass
class UnderflowError(Exception): pass
class InvalidValueError(Exception): pass


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
    elif flag == 6:
        if CATCHOVERFLOW: 
            raise OverflowError(err)
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
    elif flag == 10:
        if CATCHOVERFLOW: 
            raise OverflowError(err)
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
    """ Check the type of an object. If it thinks the type is wrong, it will tell you what the
    name of obj is (as you named it in your code) along with its type and what was expected.
    Grabbing the name doesn't work if you pass him a dictionary element like myDict['key'];
    it can only tell the name myDict. One could use type hints, but at the time of writing,
    type hints will not necessarily crash the program, which I want.

    Args:
        obj (obj)
        expectedType (type): what type do you expect? Also accepts "array", "real", and "scalar".
    """
    calling_frame = inspect.currentframe().f_back
    locals_dict = calling_frame.f_locals
    if len(locals_dict)==0:
        objName = 'OBJECT LITERAL'
    for var_name in locals_dict.keys():
        objName = var_name
        break 
    if expectedType=="array":
        if not isArrayLike(obj):
            if type(obj)==list or type(obj)==np.ndarray:
                logger.TBError('Received empty',type(obj),'for',objName,frame=3)
            else:
                logger.TBError('Expected array-like object for',objName,'but received',type(obj),frame=3)
    elif expectedType=="scalar":
        if isArrayLike(obj) or obj is None:
            logger.TBError('Expected scalar-like object for',objName,'but received',type(obj),frame=3)
    elif expectedType=="real":
        if isArrayLike(obj) or obj is None:
            logger.TBError('Expected real scalar object for',objName,'but received',type(obj),frame=3)
        elif obj.imag>0:
            logger.TBError('Expected real scalar object',objName,'has nonzero imaginary part',frame=3)
    else:
        if not isinstance(obj,expectedType):
            logger.TBError('Expected type',expectedType,'for',objName,'but received',type(obj),frame=3)


def checkDomain(obj, expectedDomain):
    """ Check that obj lies in expectedDomain.

    Args:
        obj (obj)
        expectedDomain (array-like): collection of values obj is allowed to take
    """
    calling_frame = inspect.currentframe().f_back
    locals_dict = calling_frame.f_locals
    for var_name, _ in locals_dict.items():
        objName = var_name  
    if not obj in expectedDomain:
        logger.TBError('Expected',objName,'to be one of',expectedDomain,frame=3)


def checkEqualLengths(*args):
    """ Check that all array-like objects passed have the same length. """
    length = len(envector(args[0]))
    for i in range(len(args)):
        if args[i] is not None:
            if len(envector(args[i])) != length:
                logger.info(length, len(envector(args[i])))
                logger.TBError('Array length mismatch detected on array',i,frame=3)


def checkExtension(filename,extension,ignoreExtension=False):
    """ Check the extension of a file

    Args:
        filename (str)
        extension (str)
        ignoreExtension (bool, optional): Defaults to False.
    """
    if not filename.endswith(extension) and not ignoreExtension:
        logger.TBError('Expected a',extension,'file.')
