latqcdtools.interfaces.interfaces
=============

```Python
convertTable(source, target, sourceDelimiter='', targetDelimiter=''):
'''
Convert a source table into a target table. The assumption for the source file is that
is that the only lines are table lines, i.e. there's no intervening \hline or something like that.
The table type is determined by the file extensions of source and target.

Args:
    source (str): source filename 
    target (str): target filename 
'''
```
```Python
readGPL(filename, discardTag=True, raggedWarn=True, floatT=<class 'numpy.float64'>):
'''
Load GPL files from Peter Lepage's g-2 tools as 2d array. Can also load GPL-like files, where one allows the
tag (column 0) on each line to be different. Optionally ignore tag, which is just a label. Implemented in this way
rather than using genfromtxt to allow the possibility of ragged tables. 
'''
```
```Python
readJSON(filename, ignoreExtension=False) -> dict:
'''
Load a JSON file. Returns a dict, where each key level corresponds to an organizational level of the JSON. 
'''
```
```Python
readPickle(filename):
'''
Load a Pickle file.

Args:
    filename (str)
'''
```
```Python
readWML(filename) -> list:
'''
Does its best to read a table from Wikipedia Markup Language. Returns a list of lists,
where each row corresponds to either a line of the table or a line of markup code. You
will have to do some processing by hand, since so many people edit Wikipedia and have
inconsistent styles.

Args:
    filename (str): Name of file 

Returns:
    list: list of rows and commands in markup table 
'''
```
```Python
readYAML(filename, ignoreExtension=False) -> dict:
'''
Load a YAML file. Returns a dict, where each key level corresponds to an organizational level of the YAML. 
'''
```
```Python
writeJSON(data, filename):
'''
Write dictionary to JSON file.

Args:
    data (dict)
    filename (str)
'''
```
```Python
writeYAML(data, filename):
'''
Write dictionary to YAML file.

Args:
    data (dict)
    filename (str)
'''
```
```Python
class csvTable(delimiter):
```
```Python
class genericTable(delimiter=None, pre='', post=''):
```
```Python
class latexTable():
```
```Python
class markdownTable():
```
```Python
class redmineTable():
```
