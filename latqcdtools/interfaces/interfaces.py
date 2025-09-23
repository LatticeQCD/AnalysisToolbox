# 
# interfaces.py
# 
# D. Clarke
# 
# Some common classes and functions that may be shared among multiple interfaces modules.
#

import yaml, json, pickle
import numpy as np
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.base.check import checkType, checkExtension
from latqcdtools.base.utilities import substringBetween
import latqcdtools.base.logger as logger


class HotQCD_MILC_Params(latticeParams):
    """ 
    A class to handle and check the input parameters of a lattice run using conventions common to both the
    HotQCD and MILC collaborations. 
    """

    def __repr__(self) -> str:
        return "HotQCD_MILC_Params"

    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        if self.Nf=='211' or self.Nf=='111':
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2+'m'+self.cm3
        else:
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2


def paramFromEnsLabel(ensemble,format='MILC'):
    """ 
    Given an ensemble string, get the parameters out of it. 

    Args:
        ensemble (str): ensemble label

    Returns:
        tuple: Ns, Nt, Nf, beta string, mass1 string, mass2 string, mass3 string
    """
    checkType(str,ensemble=ensemble)
    checkType(str,format=format)
    Ns, Nt, Nf, cbeta, cm1, cm2, cm3 = None, None, None, None, None, None, None
    if format=='MILC':
        NsNt = substringBetween(ensemble,'l','f') 
        if len(NsNt)==3:
            Ns=NsNt[:2]
            Nt=NsNt[-1]
        elif len(NsNt)==4:
            Ns=NsNt[:2]
            Nt=NsNt[2:]
        else:
            logger.TBRaise('I assume Ns has 2 digits and Nt has 1-2 digits.') 
        Nf    = substringBetween(ensemble,'f','b') 
        cbeta = substringBetween(ensemble,'b','m') 
        cm1   = ensemble.split('m')[1].strip()
        if len(Nf)>1:
            cm2 = ensemble.split('m')[2].strip()
        if len(Nf)>2:
            cm3 = ensemble.split('m')[3].strip()
    else:
        logger.TBRaise('Unsupported format',format)
    return int(Ns), int(Nt), Nf, cbeta, cm1, cm2, cm3 


def readGPL(filename,discardTag=True,raggedWarn=True,floatT=np.float64):
    """ 
    Load GPL files from Peter Lepage's g-2 tools as 2d array. Can also load GPL-like files, where one allows the
    tag (column 0) on each line to be different. Optionally ignore tag, which is just a label. Implemented in this way
    rather than using genfromtxt to allow the possibility of ragged tables. 
    """
    checkType(str,filename=filename)
    checkType(bool,discardTag=discardTag)
    checkType(bool,raggedWarn=raggedWarn)
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
        return np.array(data,dtype=floatT)
    else:
        return np.array(data,dtype=object)


def readPickle(filename):
    """
    Load a Pickle file.

    Args:
        filename (str)
    """
    checkType(str,filename=filename)
    try:
        with open(filename, 'rb') as file:
            data = pickle.load(file)
    except Exception as e:
        logger.warn('Is',filename,'actually a pickle file?')
        raise e 
    return data


def readYAML(filename,ignoreExtension=False) -> dict:
    """ 
    Load a YAML file. Returns a dict, where each key level corresponds to an organizational level of the YAML. 
    """
    checkType(str,filename=filename)
    if not ignoreExtension:
        checkExtension(filename,'yaml')
    with open(filename, 'r') as file:
        return yaml.safe_load(file)


def readJSON(filename,ignoreExtension=False) -> dict:
    """ 
    Load a JSON file. Returns a dict, where each key level corresponds to an organizational level of the JSON. 
    """
    checkType(str,filename=filename)
    if not ignoreExtension:
        checkExtension(filename,'json')
    with open(filename, 'r') as file:
        return json.load(file)


def readWML(filename) -> list:
    """ 
    Does its best to read a table from Wikipedia Markup Language. Returns a list of lists,
    where each row corresponds to either a line of the table or a line of markup code. You
    will have to do some processing by hand, since so many people edit Wikipedia and have
    inconsistent styles.

    Args:
        filename (str): Name of file 

    Returns:
        list: list of rows and commands in markup table 
    """
    checkType(str,filename=filename)
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
    """ 
    Write dictionary to YAML file.

    Args:
        data (dict)
        filename (str)
    """
    checkType(dict,data=data)
    checkType(str,filename=filename)
    with open(filename, 'w') as file:
        yaml.safe_dump(data, file) 


def writeJSON(data,filename):
    """ 
    Write dictionary to JSON file.

    Args:
        data (dict)
        filename (str)
    """
    checkType(dict,data=data)
    checkType(str,filename=filename)
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)


class genericTable(list):
    
    def __init__(self,delimiter=None,pre='',post=''):
        """ 
        genericTable objects are to streamline output tables in an arbitrary format. A genericTable
        is implemented as a list of lists.

        Args:
            delimiter (str)
            pre (str, optional): String to appear at beginning of every line table. Defaults to ''.
            post (str, optional): String to appear at end of every line of table. Defaults to ''.
        """
        if delimiter is None:
            logger.TBRaise("Please set a delimiter. Use delimiter='' for generic whitespace.")
        checkType(str,delimiter=delimiter)
        checkType(str,pre=pre)
        checkType(str,post=post)
        self.delimiter=delimiter
        self.pre=pre
        self.post=post

    def __repr__(self) -> str:
        return "genericTable"

    def __str__(self) -> str:
        result = '\n'
        for row in self:
            result += str(row) + '\n'
        return result 

    def append(self, item):
        checkType(list,item=item)
        super(genericTable, self).append(item)

    def empty(self):
        while len(self)>0:
            self.pop()

    def outputTable(self,filename=None):
        """ 
        Lets you output a table.

        Args:
            filename (str, optional): If set, will output to file. Otherwise output to screen. 
        """
        if filename is not None:
            checkType(str,filename=filename)
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
        """ 
        Convert a line of the table to a list. 
        """ 
        cols=line.strip().rstrip('\\').split(self.delimiter)
        row = []
        for col in cols:
            row.append(col.strip())
        self.append(row)

    def readTable(self,filename):
        """ 
        Read a table containing data only. (No headers!)

        Args:
            filename (str)
        """
        checkType(str,filename=filename)
        self.empty()
        inFile = open(filename,'r')
        for line in inFile:
            self.readLine(line)
        inFile.close()


class csvTable(genericTable):
    def __init__(self,delimiter):
        super().__init__(delimiter=delimiter)
    def __repr__(self) -> str:
        return "csvTable"


class latexTable(genericTable):
    def __init__(self):
        super().__init__(delimiter='&', pre='', post='\\\\')
    def __repr__(self) -> str:
        return "latexTable"


class markdownTable(genericTable):
    def __init__(self):
        super().__init__(delimiter='|', pre='|', post='|')
    def __repr__(self) -> str:
        return "markdownTable"
    def readLine(self,line):
        """ Convert a line of the table to a list. """ 
        cols=line.strip().split(self.delimiter)
        del cols[0]
        del cols[-1]
        row = []
        for col in cols:
            row.append(col.strip())
        self.append(row)


# I didn't realize at first that Redmine is using Markdown, and now I want to
# maintain backwards compatibility.
class redmineTable(markdownTable):
    def __repr__(self) -> str:
        return "redmineTable"


def convertTable(source,target,sourceDelimiter='',targetDelimiter=''):
    r""" 
    Convert a source table into a target table. The assumption for the source file is that
    is that the only lines are table lines, i.e. there's no intervening \hline or something like that.
    The table type is determined by the file extensions of source and target.

    Args:
        source (str): source filename 
        target (str): target filename 
    """
    checkType(str,source=source)
    checkType(str,target=target)
    checkType(str,sourceDelimiter=sourceDelimiter)
    checkType(str,targetDelimiter=targetDelimiter)
    sourceType = source.split('.')[-1]
    targetType = target.split('.')[-1]
    inFile = open(source,'r')
    if sourceType == 'tex':
        sourceTable = latexTable()
    elif sourceType == 'redmine':
        sourceTable = redmineTable()
    elif sourceType == 'md':
        sourceTable = markdownTable()
    elif sourceType == 'csv':
        sourceTable = csvTable(delimiter=sourceDelimiter)
    else: 
        logger.TBRaise('Unknown source file type',sourceType)
    if targetType == 'tex':
        targetTable = latexTable()
    elif targetType == 'redmine':
        targetTable = redmineTable()
    elif targetType == 'md':
        targetTable = markdownTable()
    elif targetType == 'csv':
        targetTable = csvTable(delimiter=targetDelimiter)
    else: 
        logger.TBRaise('Unknown target file type',targetType)
    for row in inFile:
        start = len(sourceTable.pre)
        end = len(sourceTable.post)
        if end==2: # This is for LaTeX
            end=3
        if start==end==0:
            items = row
        else:
            items = row[start:-end]
        if sourceTable.delimiter=='':
            items = items.split()
        else:
            items = items.split(sourceTable.delimiter)
        targetTable.append(items)
    targetTable.outputTable(target)
    inFile.close()
