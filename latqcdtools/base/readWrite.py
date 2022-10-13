# 
# readWrite.py
#
# H. Sandmeyer, D. Clarke
#
# Collection of methods for reading in and writing out data for different purposes. There are some methods
# specifically for reading in correlator fitting data, but also some more general methods
# that anyone can use, like read_in. 
#

import numpy as np
import latqcdtools.base.logger as logger

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


def read_in_pure(filename, col1=1, col2=2, symmetrize = False):
    xdata, data, nconfs = read_in_pure_no_numpy(filename, col1, col2, symmetrize)
    return np.array(xdata), np.array(data), nconfs


def read_in(filename, *args):
    """ General wrapper for reading in specified columns from a file. Comment character #. """
    if args == ():
        args = (1, 2, 3)
    data = [[]]*len(args)
    try:
        # To support input file streams
        ins = open(filename, "r")
        close = True
    except TypeError:
        close = False
        ins = filename
    for line in ins:
        if line.startswith("#"):
            continue
        line_elems = line.split()
        if len(line_elems) >= len(args):
            for i, col in enumerate(args):
                data[i].append(float(line_elems[col - 1]))
    if close:
        ins.close()
    return np.array(data)


def writeTable(filename,*args,**kwargs):
    """ Wrapper for np.savetxt, which would otherwise output in a way that is not very intuitive for tables.
        Additionally constructs format string to include only 8 digits after the decimal point. """
    if 'header' in kwargs:
        head = kwargs['header']
        if isinstance(head,list):
            form = '%12s'
            temp = (head[0],)
            if len(head[0]) > 12:
                logger.warn("writeTable header[0] should be kept under 12 characters.")
            for label in head[1:]:
                if len(label)>14:
                    logger.warn("writeTable header labels should be kept under 14 characters.")
                form += '  %14s'
                temp += label,
            head = form % temp
    else:
        head = ''
    data = ()
    form = ''
    for col in args:
        if isinstance(col[0],complex):
            data += (col.real,)
            data += (col.imag,)
            form += '%.8e  %8e  '
        else:
            data += (col,)
            form += '%.8e  '
    np.savetxt(filename,np.transpose(data),fmt=form,header=head)


def printClean(*args):
    data = ()
    form = ''
    for col in args:
        if isinstance(col,complex):
            data += (col.real,)
            data += (col.imag,)
            form += '(%.8e  %8e)  '
        else:
            data += (col,)
            form += '%.8e  '
    print(form % data)
