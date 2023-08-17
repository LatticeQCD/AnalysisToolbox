# 
# interfaces.py
# 
# D. Clarke
# 
# Some common classes and functions that may be shared among multiple interfaces modules.
#

import re
from sys import set_coroutine_origin_tracking_depth
import yaml
import numpy as np
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.base.check import checkType
from latqcdtools.base.utilities import substringBetween
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


def paramFrom_HotQCD_MILC(ensemble):
    checkType(ensemble,str)
    NsNt = substringBetween(ensemble,'l','f') 
    if len(NsNt)==3:
        Ns=NsNt[:2]
        Nt=NsNt[-1]
    elif len(NsNt)==4:
        Ns=NsNt[:2]
        Nt=NsNt[2:]
    else:
        logger.TBError('I do not know how to handle an ensemble name of this form.')
    Nf    = substringBetween(ensemble,'f','b') 
    cbeta = substringBetween(ensemble,'b','m') 
    cm1   = ensemble.split('m')[1].strip()
    cm2   = ensemble.split('m')[2].strip()
    return int(Ns), int(Nt), Nf, cbeta, cm1, cm2 


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
            
            
class genericTable(list):
    
    def __init__(self,delimiter=None,pre='',post=''):
        """ genericTable objects are to streamline output tables in an arbitrary format. A genericTable
        is implemented as a list of lists.

        Args:
            delimiter (str, optional)
            pre (str, optional): String to appear at beginning of every line table. Defaults to ''.
            post (str, optional): String to appear at end of every line of table. Defaults to ''.
        """
        if delimiter is None:
            logger.TBError('Please set a delimiter.')
        checkType(delimiter,str)
        checkType(pre,str)
        checkType(post,str)
        self.delimiter=delimiter
        self.pre=pre
        self.post=post

    def append(self, item):
        checkType(item, list)
        super(genericTable, self).append(item)
        
    def outputTable(self,filename=None):
        """ Lets you output a table.

        Args:
            filename (str, optional): If set, will output to file. Otherwise output to screen. 
        """
        if filename is not None:
            outFile = open(filename,'w')
        for row in self:
            line = self.pre + ' ' +  str(row[0])
            for col in range(1,len(row)):
                line += ' ' + self.delimiter + ' ' + str(row[col])
            line += ' ' + self.post
            if filename is None:
                logger.info(line)
            else:
                outFile.write(line+'\n')
        if filename is not None:
            outFile.close()


class latexTable(genericTable):
    def __init__(self):
        super().__init__(delimiter='&', pre='', post='\\\\')


class redmineTable(genericTable):
    def __init__(self):
        super().__init__(delimiter='|', pre='|', post='|')


def convertTable(source,target):
    """ Convert a source table into a target table. The assumption for the source file is that
    is that the only lines are table lines, i.e. there's no intervening \hline or something like that.
    The table type is determined by the file extensions of source and target.

    Args:
        source (str): source filename 
        target (str): target filename 
    """
    checkType(source,str)
    checkType(target,str)
    sourceType = source.split('.')[-1]
    targetType = target.split('.')[-1]
    inFile = open(source,'r')
    if sourceType == 'tex':
        sourceTable = latexTable()
    elif sourceType == 'redmine':
        sourceTable = redmineTable()
    else: 
        logger.TBError('Unknown source file type',sourceType)
    if targetType == 'tex':
        targetTable = latexTable()
    elif targetType == 'redmine':
        targetTable = redmineTable()
    else: 
        logger.TBError('Unknown target file type',targetType)
    for row in inFile:
        start = len(sourceTable.pre)
        end = len(sourceTable.post)
        if end==2: # This is for LaTeX
            end=3
        items = row[start:-end]
        items = items.split(sourceTable.delimiter)
        targetTable.append(items)
    targetTable.outputTable(target)
    inFile.close()