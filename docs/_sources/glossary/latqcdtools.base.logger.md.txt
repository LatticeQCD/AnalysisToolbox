latqcdtools.base.logger
=============

```Python
TBError(*args, frame=2):
'''
Print error message and exit with -1.

Args:
    frame (int, optional): Controls the name of the caller. Defaults to method that called TBError.
'''
```
```Python
TBFail(*args)
```
```Python
TBPass(*args)
```
```Python
TBRaise(*args, frame=2, exception=None):
'''
Add custom text to exceptions, which can be useful for debugging sometimes. 

Args:
    frame (int, optional): Controls the name of the caller. Defaults to method that called TBRaise.
    exception (Exception, optional): If None, raise generic ToolboxException. Otherwise raise exception.
'''
```
```Python
_getCallerName(frame):
'''
Gets the name of the function that calls the present function. 
'''
```
```Python
_getTimeStamp():
'''
Get HH:MM:SS 
'''
```
```Python
_log(outString)
```
```Python
createLogFile(filename='Toolbox.log'):
'''
Have output sent also to a log file filename. If this file already exists, it will get deleted. We use the
logging module because it knows how to handle multiple processes writing to the same file. 
'''
```
```Python
debug(*args, frame=2)
```
```Python
details(*args)
```
```Python
info(*args)
```
```Python
progress(*args)
```
```Python
set_log_level(level)
```
```Python
warn(*args, frame=2)
```
```Python
class ToolboxException:
```
