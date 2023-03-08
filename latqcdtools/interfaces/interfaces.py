# 
# interfaces.py
# 
# D. Clarke
# 
# Some common classes and functions that may be shared among multiple interfaces modules.
#

import yaml
import numpy as np
from latqcdtools.physics.lattice_params import latticeParams
import latqcdtools.base.logger as logger


class HotQCD_MILC_Params(latticeParams):
    """ A class to handle and check the input parameters of a lattice run using conventions common to both the
        HotQCD and MILC collaborations. """

    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        if self.Nf=='211':
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2+'m'+self.cm3
        else:
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2


def loadGPL(filename,discardTag=True):
    """ Load GPL files from Peter Lepage's g-2 tools as 2d array. Can also load GPL-like files, where one allows the
    tag (column 0) on each line to be different. Optionally ignore tag, which is just a label. Implemented in this way
    rather than using genfromtxt to allow the possibility of ragged tables. """
    gplFile = open(filename,'r')
    minIndex = 0
    data = []
    if discardTag:
        minIndex = 1
    colLengths = []
    for line in gplFile:
        parse = line.split()
        colLengths.append(len(parse))
    gplFile.close()
    gplFile = open(filename,'r')
    minLength = min(colLengths)
    maxLength = max(colLengths)
    if minLength != maxLength:
        logger.warn('Loaded ragged table. Using minLength =',minLength,'and truncating the rest.')
    for line in gplFile:
        parse = line.split()
        data.append( [ parse[i] for i in range(minIndex,minLength) ] )
    gplFile.close()
    if discardTag:
        return np.array(data,dtype=float)
    else:
        return np.array(data,dtype=object)


def loadYAML(filename):
    """ Load a YAML file. Returns a dict, where each key level corresponds to an organizational level of the YAML. """
    if not filename.endswith('yaml'):
        logger.TBError('Expected a yaml file.')
    with open(filename, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            logger.TBError('Encountered exception:',exc)