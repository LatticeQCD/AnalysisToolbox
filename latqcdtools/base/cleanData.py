# 
# cleanData.py
# 
# D. Clarke
# 
# Some methods for cleaning up numpy arrays. You might use this e.g. to restrict data to certain intervals, or to
# exclude NaN and Inf. These are essentially wrappers for array masks, which I find hard to read.
#

import numpy as np

def clipRange(array,minVal=-np.inf,maxVal=np.inf):
    """ Throw out any elements of array that lie outside the interval [minVal,maxVal].  """
    mask = np.logical_and( array[:]>=minVal, array[:]<=maxVal )
    return array[mask]


def intersectAtCol(table1,table2,col):
    """ Return only those rows of left and right that have identical elements in column col. """
    table1       = np.array(table1)
    table2       = np.array(table2)
    mask1using2  = np.isin( table1[col,:], table2[col,:] )
    mask2using1  = np.isin( table2[col,:], table1[col,:] )
    return table1[:,mask1using2], table2[:,mask2using1]
