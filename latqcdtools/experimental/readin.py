# 
# readin.py                                                               
# 
# Collection of methods for reading in data for different purposes. There are some methods
# specifically for reading in correlator fitting data, but also some more general methods
# that anyone can use, like read_in. 
#
 
import numpy as np


def read_in_pure_no_numpy(filename, col1=1, col2=2, symmetrize = False):
    try:
        #To support input file streams
        ins = open(filename, "r")
        close = True
    except TypeError:
        ins = filename
        close = False
    data_dict = {}
    lineelems = []
    line_numb = 0
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


def read_in_pure_trans(filename, col1=1, col2=2, symmetrize = False):
    xdata, data, nconfs = read_in_pure(filename, col1, col2, symmetrize)
    return xdata, np.array(data).transpose(), len(data[0])


''' General wrapper for reading in specified columns from a file. Comment character #. '''
def read_in(filename, *args):
    if args == ():
        args = (1, 2, 3)
    data = [[] for i in range(len(args))]
    try:
        #To support input file streams
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


def read_in_fitmass(filename):

    with open(filename) as f:
        header = f.readline()
    
    cols = header.split()

    ns_cols = []
    Ans_cols = []
    for ns in range(100):
        if header.find("m_no_" + str(ns)) == -1:
            break
        else:
            ns_cols.append(cols.index("m_no_" + str(ns)))
            Ans_cols.append(cols.index("A_no_" + str(ns)))

    nstates = ns

    nso_cols = []
    Anso_cols = []
    for nso in range(100):
        if header.find("m_osc_" + str(nso)) == -1:
            break
        else:
            nso_cols.append(cols.index("m_osc_" + str(nso)))
            Anso_cols.append(cols.index("A_osc_" + str(nso)))

    nstates_osc = nso

    nparams = 2 * (nstates + nstates_osc)

    data = np.loadtxt(filename)

    aicc_col = cols.index("AICc")
    chi2_col = cols.index("chi^2/d.o.f.")

    #When there is only one line, numpy.loadtxt, will giv a 1D array
    if len(data.shape) == 1:
        data = [data]

    ranges = []
    res = []
    res_err = []
    aicc = []
    chi2_dof = []

    for line in data:
        ranges.append((line[0], line[-1]))
        aicc.append(line[aicc_col])
        chi2_dof.append(line[chi2_col])
        res.append(np.ones(nparams))
        res_err.append(np.ones(nparams))
        
        for ns in range(nstates):
            res[-1][2 * ns + 1] = line[ns_cols[ns]]
            res_err[-1][2 * ns + 1] = line[ns_cols[ns] + 1]
        for nso in range(nstates_osc):
            res[-1][2 * nstates + 2 * nso + 1] = line[nso_cols[nso]]
            res_err[-1][2 * nstates + 2 * nso + 1] = line[nso_cols[nso] + 1]

        for ns in range(nstates):
            res[-1][2 * ns] = line[Ans_cols[ns]]
            res_err[-1][2 * ns] = line[Ans_cols[ns] + 1]
        for nso in range(nstates_osc):
            res[-1][2 * nstates + 2 * nso] = line[Anso_cols[nso]]
            res_err[-1][2 * nstates + 2 * nso] = line[Anso_cols[nso] + 1]

    return (np.array(ranges, dtype = int), np.array(res), np.array(res_err), np.array(aicc),
                np.array(chi2_dof), np.array(nstates), np.array(nstates_osc))
