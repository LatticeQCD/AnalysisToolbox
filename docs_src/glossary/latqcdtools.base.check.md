latqcdtools.base.check
=============

`checkDomain(obj, expectedDomain)`
 
    Check that obj lies in expectedDomain.

    Args:
        obj (obj)
        expectedDomain (array-like): collection of values obj is allowed to take
    
`checkEqualLengths(*args)`
 
    Check that all array-like objects passed have the same length. 
    
`checkExtension(filename, extension, ignoreExtension=False)`
 
    Check the extension of a file

    Args:
        filename (str)
        extension (str)
        ignoreExtension (bool, optional): Defaults to False.
    
`checkType(expectedType, **kwargs)`
 
    Check the type of an object. If it thinks the type is wrong, it will tell you what the
    name of obj is (as you named it in your code) along with its type and what was expected.
    Grabbing the name doesn't work if you pass him a dictionary element like myDict['key'];
    it can only tell the name myDict. One could use type hints, but at the time of writing,
    type hints will not necessarily crash the program, which I want.

    Args:
        obj (obj)
        expectedType (type): what type do you expect? Also accepts "array", "real", "int", and "scalar".
    
`err_handler(err, flag)`
 
    This method lets us control in detail how different types of errors are treated. 
    
`ignoreDivideByZero()`
 
    Turn off zero division crashes. 
    
`ignoreInvalidValue()`
 
    Turn off invalid value crashes. 
    
`ignoreOverflow()`
 
    Turn off overflow crashes. 
    
`ignoreUnderflow()`
 
    Turn off underflow crashes. 
    
`DivideByZeroError
`


`InvalidValueError
`


`UnderflowError
`


