latqcdtools.base.logger
=============

```python
TBError(*args, frame=2):
"""
Print error message and exit with -1.

Args:
    frame (int, optional): Controls the name of the caller. Defaults to method that called TBError.
"""
```
```python
TBFail(*args)
```
```python
TBPass(*args)
```
```python
TBRaise(*args, frame=2, exception=None):
"""
Add custom text to exceptions, which can be useful for debugging sometimes. 

Args:
    frame (int, optional): Controls the name of the caller. Defaults to method that called TBRaise.
    exception (Exception, optional): If None, raise generic ToolboxException. Otherwise raise exception.
"""
```
```python
_getCallerName(frame):
"""
Gets the name of the function that calls the present function. 
"""
```
```python
_getTimeStamp():
"""
Get HH:MM:SS 
"""
```
```python
_log(outString)
```
```python
createLogFile(filename='Toolbox.log'):
"""
Have output sent also to a log file filename. If this file already exists, it will get deleted. We use the
logging module because it knows how to handle multiple processes writing to the same file. 
"""
```
```python
debug(*args, frame=2)
```
```python
details(*args)
```
```python
info(*args)
```
```python
progress(*args)
```
```python
set_log_level(level)
```
```python
warn(*args, frame=2)
```
```python
class ToolboxException:
```
