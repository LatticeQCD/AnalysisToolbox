latqcdtools.interfaces.interfaces
=============

```python
convertTable(source, target, sourceDelimiter='', targetDelimiter=''):
"""
Convert a source table into a target table. The assumption for the source file is that
is that the only lines are table lines, i.e. there's no intervening \hline or something like that.
The table type is determined by the file extensions of source and target.

Args:
    source (str): source filename 
    target (str): target filename 
"""
```
```python
readGPL(filename, discardTag=True, raggedWarn=True, floatT=numpy.float64):
"""
Load GPL files from Peter Lepage's g-2 tools as 2d array. Can also load GPL-like files, where one allows the
tag (column 0) on each line to be different. Optionally ignore tag, which is just a label. Implemented in this way
rather than using genfromtxt to allow the possibility of ragged tables. 
"""
```
```python
readJSON(filename, ignoreExtension=False) -> dict:
"""
Load a JSON file. Returns a dict, where each key level corresponds to an organizational level of the JSON. 
"""
```
```python
readPickle(filename):
"""
Load a Pickle file.

Args:
    filename (str)
"""
```
```python
readWML(filename) -> list:
"""
Does its best to read a table from Wikipedia Markup Language. Returns a list of lists,
where each row corresponds to either a line of the table or a line of markup code. You
will have to do some processing by hand, since so many people edit Wikipedia and have
inconsistent styles.

Args:
    filename (str): Name of file 

Returns:
    list: list of rows and commands in markup table 
"""
```
```python
readXML(filename, ignoreExtension=False) -> dict:
"""
Load a XML file. Returns a dict, where each key level corresponds to an organizational level of the XML. 
"""
```
```python
readYAML(filename, ignoreExtension=False) -> dict:
"""
Load a YAML file. Returns a dict, where each key level corresponds to an organizational level of the YAML. 
"""
```
```python
writeJSON(data, filename):
"""
Write dictionary to JSON file.

Args:
    data (dict)
    filename (str)
"""
```
```python
writeYAML(data, filename):
"""
Write dictionary to YAML file.

Args:
    data (dict)
    filename (str)
"""
```
```python
class csvTable(delimiter):
```
```python
class genericTable(delimiter=None, pre='', post=''):
```
```python
class latexTable():
```
```python
class markdownTable():
```
```python
class redmineTable():
```
