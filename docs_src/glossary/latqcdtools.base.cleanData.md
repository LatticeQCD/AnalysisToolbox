latqcdtools.base.cleanData
=============

`clipRange(array, col=None, minVal=-inf, maxVal=inf) -> numpy.ndarray`
 
    Throw out any elements of array that lie outside the interval (minVal,maxVal). Note this
    renders arrays finite. 
    
`deleteCol(array, col) -> numpy.ndarray`


`deleteRow(array, row) -> numpy.ndarray`


`excludeAtCol(table, col, atVal) -> numpy.ndarray`
 
    Return everything except those rows of table where col has exactly the value atVal. 
    
`intersectAtCol(table1, table2, col) -> numpy.ndarray`
 
    Return only those rows of table1 and table2 that have identical elements in column col. 
    
`restrictAtCol(table, col, atVal) -> numpy.ndarray`
 
    Return only those rows of table where col has exactly the value atVal. 
    
`spliceAtCol(table1, table2, col, atVal) -> numpy.ndarray`
 
    Assuming two tables table1 and table2 have common values in column col, create a new
    table, where table1 has corresponding entries less than atVal in col, and table 2
    has corresponding entries greater than atVal. 
    
