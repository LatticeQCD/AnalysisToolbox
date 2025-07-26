latqcdtools.base.cleanData
=============

```Python
clipRange(array, col=None, minVal=-inf, maxVal=inf) -> numpy.ndarray:
'''
Throw out any elements of array that lie outside the interval (minVal,maxVal). Note this
renders arrays finite. 
'''
```
```Python
deleteCol(array, col) -> numpy.ndarray:
'''
Remove a column of a 2d np.ndarray.

Args:
    array (np.ndarray)
    col (int): Remove this column 

Returns:
    np.ndarray: array with row removed
'''
```
```Python
deleteRow(array, row) -> numpy.ndarray:
'''
Remove a row of a 2d np.ndarray.

Args:
    array (np.ndarray)
    row (int): Remove this row 

Returns:
    np.ndarray: array with row removed
'''
```
```Python
excludeAtCol(table, col=None, atVal=inf) -> numpy.ndarray:
'''
Return everything except those rows of table where col has exactly the value atVal. 
'''
```
```Python
intersectAtCol(table1, table2, col):
'''
Return only those rows of table1 and table2 that have identical elements in column col. 
'''
```
```Python
restrictAtCol(table, col, atVal, rtol=None, atol=None) -> numpy.ndarray:
'''
Return only those rows of table where col has exactly the value atVal. 
'''
```
```Python
spliceAtCol(table1, table2, col, atVal) -> numpy.ndarray:
'''
Assuming two tables table1 and table2 have common values in column col, create a new
table, where table1 has corresponding entries less than atVal in col, and table 2
has corresponding entries greater than atVal. 
'''
```
