latqcdtools.base.check
=============

```Python
checkDomain(obj, expectedDomain):
'''
Check that obj lies in expectedDomain.

Args:
    obj (obj)
    expectedDomain (array-like): collection of values obj is allowed to take
'''
```
```Python
checkEqualLengths(*args):
'''
Check that all array-like objects passed have the same length. 
'''
```
```Python
checkExtension(filename, extension):
'''
Check the extension of a file

Args:
    filename (str)
    extension (str)
'''
```
```Python
checkType(expectedType, **kwargs):
'''
Check the type of an object. If it thinks the type is wrong, it will tell you what the
name of obj is (as you named it in your code) along with its type and what was expected.
Grabbing the name doesn't work if you pass him a dictionary element like myDict['key'];
it can only tell the name myDict. One could use type hints, but at the time of writing,
type hints will not necessarily crash the program, which I want.

Args:
    obj (obj)
    expectedType (type): what type do you expect? Also accepts "array", "real", "int", and "scalar".
'''
```
```Python
err_handler(err, flag):
'''
This method lets us control in detail how different types of errors are treated. 
'''
```
```Python
ignoreDivideByZero():
'''
Turn off zero division crashes. 
'''
```
```Python
ignoreInvalidValue():
'''
Turn off invalid value crashes. 
'''
```
```Python
ignoreOverflow():
'''
Turn off overflow crashes. 
'''
```
```Python
ignoreUnderflow():
'''
Turn off underflow crashes. 
'''
```
```Python
class DivideByZeroError:
```
```Python
class InvalidValueError:
```
```Python
class UnderflowError:
```
