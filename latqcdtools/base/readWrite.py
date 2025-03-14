# 
# readWrite.py
#
# D. Clarke
#
# Methods for convenient reading and writing. 
#


import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.base.cleanData import clipRange, excludeAtCol 
from latqcdtools.base.utilities import createFilePath


def readTable(filename,unpack=True,col=None,minVal=-np.inf,maxVal=np.inf,excludeAtVal=None,**kwargs) -> np.ndarray:
    """ 
    Wrapper for np.loadtxt. It unpacks by default to prevent transposition errors, and also optionally
    allows the user to restrict the table based on the range of one of the columns.

    Args:
        filename (str):
            Name of the file to open. 
        unpack (bool, optional): 
            If False, reads the table in a confusing, useless way. Defaults to True.
        col (int, optional): 
            Use this column to restrict the range of the rest of the table. Defaults to None.
        minVal (float, optional): 
            Minimum value for above restriction. Defaults to -np.inf.
        maxVal (float, optional):
            Maximum value for above restriction. Defaults to np.inf.
        excludeAtVal (float, optional):
            Throw out all rows for which column col has value excludeAtVal. Exclusion is
            carried out after range restriction. 
            
    Returns:
        np.array: Data table. 
    """
    checkType(str,filename=filename)
    try: 
        data = np.loadtxt(filename,unpack=unpack,**kwargs)
    except Exception as e:
        raise e
    if col is not None:
        data = clipRange(data,col=col,minVal=minVal,maxVal=maxVal)
        if excludeAtVal is not None:
            data = excludeAtCol(data,col=col,atVal=excludeAtVal)
    else:
        if minVal!=-np.inf:
            logger.TBRaise('Set col=None with minVal')
        if maxVal!=np.inf:
            logger.TBRaise('Set col=None with maxVal')
        if excludeAtVal is not None:
            logger.TBRaise('Set col=None with excludeAtVal')
    return data


def writeTable(filename,*args,**kwargs):
    """ 
    Wrapper for np.savetxt. The idea is that you can use it like this:
    
    writeTable('file.txt',col1,col2,header=['header1','header2])
    
    This works for an arbitrary number of 1-d columns col. It seems much more intuitive to me 
    that you pass columns as arguments than whatever np.savetxt is doing.

    Args:
        filename (str): output file name
    """
    npkwargs=kwargs
    if len(args)==0:
        logger.TBRaise('No data passed to writeTable.')
    if 'header' in kwargs:
        head = kwargs['header']
        del npkwargs['header']
        if isinstance(head,list):
            form = '%15s'
            temp = (head[0],)
            if len(head[0]) > 12:
                logger.warn("writeTable header[0] should be kept under 12 characters.")
            for label in head[1:]:
                if len(label)>15:
                    logger.warn("writeTable header labels should be kept under 14 characters.")
                form += '  %15s'
                temp += label,
            head = form % temp
    else:
        head = ''
    data = ()
    dtypes = []
    form = ''
    colno = 0
    ndat = len(args[0])
    for col in args:
        col_arr = np.array(col)
        if len(col_arr) != ndat:
            logger.TBRaise('Expected length',ndat,'for col',colno,'but found',len(col_arr))
        if isinstance(col_arr[0],complex):
            data += (col_arr.real,)
            data += (col_arr.imag,)
            form += '  %15.8e  %15.8e'
            dtypes.append( (_lab(colno), float) )
            dtypes.append( (_lab(colno+1), float) )
            colno += 2
        elif isinstance(col_arr[0],str):
            data += (col_arr,)
            form += '  %15s'
            dtypes.append( (_lab(colno), 'U15' ) ) # 15 characters
            colno += 1
        else:
            data += (col_arr,)
            form += '  %15.8e'
            dtypes.append( (_lab(colno), float) )
            colno += 1
    ab = np.zeros(data[0].size, dtype=dtypes)
    for i in range(colno):
        ab[_lab(i)] = data[i]
    createFilePath(filename)
    np.savetxt(filename, ab, fmt=form, header=head, **kwargs)


def _lab(num) -> str:
    """ 
    Create a short string label for each column of a data table. Needed for writeTable. 
    """
    return 'var' + str(num)
