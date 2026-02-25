latqcdtools.base.fileSystem
=============

```Python
createFilePath(target):
'''
Create the directory path for a file if it isn't there already. 

Args:
    target (str)
'''
```
```Python
deleteLine(target, line_number):
'''
Delete line line_number from file target, indexed from 1.

Args:
    target (str)
    line_number (int)
'''
```
```Python
getFileTimeStamp(target) -> str:
'''
Get the time stamp (when it was last modified) of a regular file.

Args:
    target (str)

Returns:
    str: time stamp in format 2025-09-23 14:56:27
'''
```
```Python
ls(target) -> list:
'''
Get list of files in file path. Similar to ls in bash.

Args:
    target (str)

Returns:
    list: List of files in target 
'''
```
```Python
rm(target):
'''
Delete target regular file or folder. Equivalent to rm -rf in bash.

Args:
    target (str)
'''
```
