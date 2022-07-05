# 
# cleanData.py
# 
# D. Clarke
# 
# Some methods for cleaning up numpy arrays. You might use this e.g. to restrict data to certain intervals, or to
# exclude NaN and Inf. These are essentially wrappers for array masks, which I find hard to read.
#

import numpy as np


def clipRange(array, col = None, minVal=-np.inf, maxVal=np.inf):
    """ Throw out any elements of array that lie outside the interval [minVal,maxVal].  """
    if col is None:
        mask = np.logical_and( array[:]>=minVal, array[:]<=maxVal )
        return array[mask]
    else:
        mask = np.logical_and( array[col,:]>=minVal, array[col,:]<=maxVal )
        return array[:,mask]


def intersectAtCol(table1, table2, col):
    """ Return only those rows of table1 and table2 that have identical elements in column col. """
    table1       = np.array(table1)
    table2       = np.array(table2)
    mask1using2  = np.isin( table1[col,:], table2[col,:] )
    mask2using1  = np.isin( table2[col,:], table1[col,:] )
    return table1[:,mask1using2], table2[:,mask2using1]


def spliceAtCol(table1, table2, col, atVal):
    """ Assuming two tables table1 and table2 have common values in column col, create a new
        table, where table1 has corresponding entries less than atVal in col, and table 2
        has corresponding entries greater than atVal. """
    table1 = np.array(table1)
    table2 = np.array(table2)
    mask1  = table1[col]<=atVal
    mask2  = table2[col]>atVal
    result = np.column_stack( (table1[:,mask1], table2[:,mask2]) )
    return result


def restrictAtCol(table, col, atVal):
    """ Return only those rows of table where col has exactly the value atVal. """
    table = np.array(table)
    mask  = np.equal(table[col,:],atVal)
    return table[:,mask]


def excludeAtCol(table, col, atVal):
    """ Return everything except those rows of table where col has exactly the value atVal. """
    table = np.array(table)
    mask  = np.not_equal(table[col,:],atVal)
    return table[:,mask]
