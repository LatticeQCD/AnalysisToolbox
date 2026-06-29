latqcdtools.base.utilities
=============

```python
_alphanum_key(key):
"""
Splits the string `key` at any point where there's a change from digit to non-digit or vice versa.
Will accurately handle integers and real numbers.
"""
```
```python
_convert(text)
```
```python
_getPrefix(byteString)
```
```python
appendToDocstring(string=None, args=None, returns=None)
```
```python
byteConvert(x, b1, b2):
"""
Convert between bytes given scientific prefixes.

Args:
    x (float): Bytes in original units. 
    b1 (str): Original units. 
    b2 (str): Target units.

Returns:
    float: Bytes in target units. 
"""
```
```python
cleanOutput(*args, label=None, sspace=20) -> str:
"""
This method takes a bunch of args and formats them automatically for output. The idea is
that you can use this method to ensure that columns are well lined up.

Args:
    *args: The numbers you want to output, separated by commas. 
    label (str, optional): Put label to the left of your output. Defaults to None.

Returns:
    str: formatted output string 
"""
```
```python
comesBefore(date1, date2, format='%Y/%m/%d %H:%M:%S', beforeOrEqual=False) -> bool:
"""
Check whether date1 comes before date2.

Args:
    date1 (str)
    date2 (str)
    format (str): format for date strings. Defaults to "%Y/%m/%d %H:%M:%S"
    beforeOrEqual (bool): also return True if date1 == date2
    
Returns:
    bool: date1 < date2 
"""
```
```python
elapsedSeconds(date1, date2, format='%Y/%m/%d %H:%M:%S') -> float:
"""
Compute elapsed time in seconds between date1 and date2. 

Args:
    date1 (str)
    date2 (str)
    format (str): format for date strings. Defaults to "%Y/%m/%d %H:%M:%S"
    
Returns:
    float: elapsed time in seconds 
"""
```
```python
envector(*args):
"""
Change obj to a numpy array if it's a scalar. Sometimes required when, e.g., using np.vectorize. 
"""
```
```python
find_nearest_idx(array, value) -> int:
"""
Find the index of the element of array nearest to value. 
"""
```
```python
getArgs(parser):
"""
Get arguments from the ArgumentParser. Complain if you don't get exactly the correct arguments. 
"""
```
```python
isArrayLike(obj) -> bool:
"""
Figure out whether obj is indexable.

Args:
    obj (python object)

Returns:
    bool: True if there is at least one index, false otherwise. 
"""
```
```python
isComplexType(obj) -> bool
```
```python
isFloatType(obj) -> bool
```
```python
isHigherDimensional(obj) -> bool:
"""
Figure out whether obj has at least two indices.

Args:
    obj (array-like)

Returns:
    bool: True if there are at least two indices, false otherwise. 
"""
```
```python
isIntType(obj) -> bool
```
```python
isReal(obj) -> bool
```
```python
isScalar(obj) -> bool
```
```python
naturalSort(l) -> list:
"""
Sort list of strings so that, e.g. '10' comes after '9' rather than before it.
"""
```
```python
printArg(message, param):
"""
Some arguments are None by default, and you only want to print them if they are set. 
"""
```
```python
printClean(*args, label=None, sspace=20):
"""
Wrapper for cleanOutput that prints to screen.

Args:
    *args: The numbers you want to output, separated by commas. 
    label (str, optional): Put label to the left of your output. Defaults to None.
"""
```
```python
printDict(dic, level=0):
"""
Prints key, value pairs line by line. 
"""
```
```python
shell(*args):
"""
Carry out the passed arguments args in the shell. Can be passed as a single
string or as a list. Captures and returns output of shell command. E.g.
    shell('ls -lah')
"""
```
```python
shellVerbose(*args):
"""
Same as shell, but instead of capturing output, print it to screen. 
"""
```
```python
substringBetween(string, a, b) -> str:
"""
Find the substring of string between a and b. If a==b, it looks between the
first and second occurences of a. 

Args:
    string (str)
    a (str): starting delimiter 
    b (str): ending delimiter

Returns:
    str: substring
"""
```
```python
toNumpy(*args, **kwargs)
```
```python
unvector(obj):
"""
Remove outermost brackets of array-like object with single element, if possible. This is needed
because sometimes different numpy methods give inconsistent outputs, like turning a scalar
into a zero-dimensional array, a 1-dimensional array, or just the scalar itself.

Args:
    obj (python object)

Returns:
    obj, obj[0], or obj.item() depending on obj

"""
```
```python
class timer():
"""
A class to facilitate doing rudimentary timings in the Toolbox. 
"""
```
