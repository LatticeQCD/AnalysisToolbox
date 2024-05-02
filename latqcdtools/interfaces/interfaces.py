# 
# interfaces.py
# 
# D. Clarke
# 
# Some common classes and functions that may be shared among multiple interfaces modules.
#

import yaml, json
import numpy as np
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.base.check import checkType, checkExtension
from latqcdtools.base.utilities import substringBetween
import latqcdtools.base.logger as logger


class HotQCD_MILC_Params(latticeParams):
    """ A class to handle and check the input parameters of a lattice run using conventions common to both the
        HotQCD and MILC collaborations. """

    def __repr__(self) -> str:
        return "HotQCD_MILC_Params"

    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        if self.Nf=='211':
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2+'m'+self.cm3
        else:
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2


def paramFrom_HotQCD_MILC(ensemble):
    """ Given an ensemble string of the form used by HotQCD and MILC, get all the parameters.

    Args:
        ensemble (str): ensemble label of the form l3216f3b6050m00394m1064

    Returns:
        tuple: Ns, Nt, Nf, beta string, mass1 string, mass2 string
    """
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


def readGPL(filename,discardTag=True,raggedWarn=True):
    """ Load GPL files from Peter Lepage's g-2 tools as 2d array. Can also load GPL-like files, where one allows the
    tag (column 0) on each line to be different. Optionally ignore tag, which is just a label. Implemented in this way
    rather than using genfromtxt to allow the possibility of ragged tables. """
    checkType(filename,str)
    checkType(discardTag,bool)
    checkType(raggedWarn,bool)
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
    if (minLength != maxLength) and raggedWarn:
        logger.warn('Loaded ragged table. Using minLength =',minLength,'and truncating the rest.')
    for line in gplFile:
        parse = line.split()
        data.append( [ parse[i] for i in range(minIndex,minLength) ] )
    gplFile.close()
    if discardTag:
        return np.array(data,dtype=float)
    else:
        return np.array(data,dtype=object)


def readYAML(filename,ignoreExtension=False) -> dict:
    """ Load a YAML file. Returns a dict, where each key level corresponds to an organizational level of the YAML. """
    checkType(filename,str)
    checkExtension(filename,'yaml',ignoreExtension)
    with open(filename, 'r') as file:
        try:
            return yaml.safe_load(file)
        except yaml.YAMLError as e:
            logger.TBError('Encountered exception:',e)


def readJSON(filename,ignoreExtension=False) -> dict:
    """ Load a JSON file. Returns a dict, where each key level corresponds to an organizational level of the JSON. """
    checkType(filename,str)
    checkExtension(filename,'json',ignoreExtension)
    with open(filename, 'r') as file:
        return json.load(file)


def readWML(filename) -> list:
    """ Does its best to read a table from Wikipedia Markup Language. Returns a list of lists,
    where each row corresponds to either a line of the table or a line of markup code. You
    will have to do some processing by hand, since so many people edit Wikipedia and have
    inconsistent styles.

    Args:
        filename (str): Name of file 

    Returns:
        list: list of rows and commands in markup table 
    """
    checkType(filename,str)
    wmlFile = open(filename,'r')
    data = []
    row = []
    for line in wmlFile:
        col=line.split('|')
        cleanCol = []
        for item in col:
            cleanCol.append(item.strip())
        for item in cleanCol:
            if item=='':
                continue
            if item=='-':
                cleanRow = []
                compressing = False
                cleanItem=''
                for jtem in row:
                    cleanItem += jtem
                    if jtem.startswith('{{') or jtem.startswith('[['):
                        compressing = True
                    if jtem.endswith('}}') or jtem.endswith(']]') or jtem.endswith('}}*'):
                        compressing = False
                    if not compressing:
                        cleanRow.append(cleanItem)
                        cleanItem='' 
                if len(cleanRow)>0:
                    data.append(cleanRow)
                    row=[]
            else:
                row.append(item)
    return data


def writeYAML(data,filename):
    """ Write dictionary to YAML file.

    Args:
        data (dict)
        filename (str)
    """
    checkType(data,dict)
    checkType(filename,str)
    with open(filename, 'w') as file:
        try:
            yaml.safe_dump(data, file) 
        except yaml.YAMLError as e:
            logger.TBError('Encountered exception:',e)


def writeJSON(data,filename):
    """ Write dictionary to JSON file.

    Args:
        data (dict)
        filename (str)
    """
    checkType(data,dict)
    checkType(filename,str)
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)


class genericTable(list):
    
    def __init__(self,delimiter=None,pre='',post=''):
        """ genericTable objects are to streamline output tables in an arbitrary format. A genericTable
        is implemented as a list of lists.

        Args:
            delimiter (str)
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

    def __repr__(self) -> str:
        return "genericTable"

    def append(self, item):
        checkType(item, list)
        super(genericTable, self).append(item)

    def empty(self):
        while len(self)>0:
            self.pop()

    def outputTable(self,filename=None):
        """ Lets you output a table.

        Args:
            filename (str, optional): If set, will output to file. Otherwise output to screen. 
        """
        if filename is not None:
            checkType(filename,str)
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

    def readLine(self,line):
        """ Convert a line of the table to a list. """ 
        cols=line.strip().rstrip('\\').split(self.delimiter)
        row = []
        for col in cols:
            row.append(col.strip())
        self.append(row)

    def readTable(self,filename):
        """ Read a table containing data only. (No headers!)

        Args:
            filename (str)
        """
        checkType(filename,str)
        self.empty()
        inFile = open(filename,'r')
        for line in inFile:
            self.readLine(line)
        inFile.close()


class latexTable(genericTable):
    def __init__(self):
        super().__init__(delimiter='&', pre='', post='\\\\')
    def __repr__(self) -> str:
        return "latexTable"


class redmineTable(genericTable):
    def __init__(self):
        super().__init__(delimiter='|', pre='|', post='|')
    def __repr__(self) -> str:
        return "redmineTable"
    def readLine(self,line):
        """ Convert a line of the table to a list. """ 
        cols=line.strip().split(self.delimiter)
        del cols[0]
        del cols[-1]
        row = []
        for col in cols:
            row.append(col.strip())
        self.append(row)


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