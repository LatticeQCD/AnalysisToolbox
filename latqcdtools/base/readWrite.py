# 
# readWrite.py
#
# H. Sandmeyer, D. Clarke
#
# Methods for convenient reading and writing, some tailored to specific contexts like correlator measurements.
#


import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.base.cleanData import clipRange


def readTable(filename,unpack=True,col=None,minVal=-np.inf,maxVal=np.inf,**kwargs):
    """ Wrapper for np.loadtxt. It unpacks by default to prevent transposition errors, and also optionally
    allows the user to restrict the table based on the range of one of the columns.

    Args:
        filename (str):
            Name of the file to open. 
        unpack (bool, optional): 
            If False, reads the table in a confusing, useless way. Defaults to True.
        col (int, optional): 
            Use this column to restrict the range of the rest of the table. Defaults to None.
        minVal (float, optional): 
            Minimum value for above restriction. Defaults to None.
        maxVal (float, optional):
            Maximum value for above restriction. Defaults to None.

    Returns:
        np.array: Data table. 
    """
    checkType(filename,str)
    try: 
        data = np.loadtxt(filename,unpack=unpack,**kwargs)
    except Exception as e:
        logger.TBError('Encountered exception:',e)
    if col is not None:
        data = clipRange(data,col=col,minVal=minVal,maxVal=maxVal)
    return data


def writeTable(filename,*args,**kwargs):
    """ Wrapper for np.savetxt, which is in the author's opinion the single worst piece of trash in the entire numpy
    package. Looking at how much effort it took me to tame it, I'm not sure if it was ever worth using to being with.

    Parameters
    ----------
    filename : str
        Name of output file.
    args :
        Put in all the 1-d numpy arrays you would like to write out.
    """
    if 'header' in kwargs:
        head = kwargs['header']
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
    for col in args:
        col_arr = np.array(col)
        if isinstance(col_arr[0],complex):
            data += (col_arr.real,)
            data += (col_arr.imag,)
            form += '  %15.8e  %15.8e'
            dtypes.append( (lab(colno), float) )
            dtypes.append( (lab(colno+1), float) )
            colno += 2
        elif isinstance(col_arr[0],str):
            data += (col_arr,)
            form += '  %12s'
            dtypes.append( (lab(colno), 'U12' ) ) # 12 characters
            colno += 1
        else:
            data += (col_arr,)
            form += '  %15.8e'
            dtypes.append( (lab(colno), float) )
            colno += 1
    ab = np.zeros(data[0].size, dtype=dtypes)
    for i in range(colno):
        ab[lab(i)] = data[i]
    np.savetxt(filename, ab, fmt=form, header=head)


def lab(num):
    """ Create a short string label for each column of a data table. Needed for writeTable. """
    return 'var' + str(num)


def read_in_pure_no_numpy(filename, col1=1, col2=2, symmetrize = False):
    try:
        # To support input file streams
        ins = open(filename, "r")
        close = True
    except TypeError:
        ins = filename
        close = False
    data_dict = {}
    for line in ins:
        if line.startswith('#') or len(line) < 2:
            continue
        lineelems = line.strip().split()
        try:
            Nt = int(lineelems[col1 - 1])
        except ValueError:
            Nt = float(lineelems[col1 - 1])
        corr = float(lineelems[col2 - 1])
        if Nt not in data_dict:
            data_dict[Nt] = []
        data_dict[Nt].append(corr)

    xdata = list(sorted(data_dict))
    data = [ data_dict[key] for key in sorted(data_dict.keys()) ]
    Nt = len(data)
    if symmetrize:
        if max(xdata) != Nt - 1:
            raise ValueError("The number of x values does not correspond to the largest of its values")
        if Nt % 2 != 0:
            raise ValueError("Nt must be even!")

        for i in range(len(data[0])):
            for nt in range(1, int(len(data)/2)):
                data[nt][i] = data[Nt - nt][i] = (data[nt][i] + data[Nt - nt][i]) / 2

    if close:
        ins.close()
    return xdata, data, len(data[0])


def readCorrelatorTable(filename, col1=1, col2=2, symmetrize = False):
    xdata, data, nconfs = read_in_pure_no_numpy(filename, col1, col2, symmetrize)
    return np.array(xdata), np.array(data), nconfs