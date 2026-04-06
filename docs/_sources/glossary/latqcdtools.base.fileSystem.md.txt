latqcdtools.base.fileSystem
=============

```Python
cd(target):
'''
Change to target directory. Equivalent to cd in Bash.

Args:
    target (str)
'''
```
```Python
cp(source, target):
'''
Copy source to target. Creates target directory path if needed. Similar
to cp -r in Bash.

Args:
    source (str)
    target (str)
'''
```
```Python
createFilePath(filePath):
'''
Create the directory path for a file if it isn't there already. 

Args:
    filePath (str)
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
getFileSize(target) -> int:
'''
Get size of regular file.

Args:
    target (str)

Returns:
    int: size of file in bytes 
'''
```
```Python
getFileTimeStamp(target, form='human', zone=None) -> str:
'''
Get the time stamp (when it was last modified) of a regular file.

Args:
    target (str)
    form (str)
    zone (str): list of time zones here: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
    
Returns:
    str: time stamp in format 2025-09-23 14:56:27
'''
```
```Python
getNumberLines(target) -> int:
'''
Get number of lines in human-readable file.

Args:
    target (str)

Returns:
    int: number of lines 
'''
```
```Python
ls(target) -> list:
'''
Get list of files in file path. Similar to ls in Bash.

Args:
    target (str)

Returns:
    list: List of files in target 
'''
```
```Python
rm(target):
'''
Delete target regular file or folder. Equivalent to rm -rf in Bash.

Args:
    target (str)
'''
```
