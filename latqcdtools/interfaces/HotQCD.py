# 
# HotQCD.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the HotQCD collaboration.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.math.math import rel_check
from latqcdtools.interfaces.interfaces import HotQCD_MILC_Params
from latqcdtools.base.logger import ToolboxException


class HotQCDParams(HotQCD_MILC_Params):
    """
    A class to handle and check the input parameters of a lattice run, especially for HotQCD.
    """
    def __repr__(self) -> str:
        return "HotQCDParams"
    def getcGradFlowPureGauge(self) -> str:
        return 's'+str(self.Ns).zfill(3)+'t'+str(self.Nt).zfill(2)+'_b0'+self.cbeta+'00'


def makeConfTag(conf,stream) -> str:
    """ 
    This takes a configuration number conf and stream label stream to make a tag labelling a configuration.
    Implementing this as a function makes sure everyone using the Toolbox as the same convention. 
    """
    return f'{stream}:{conf}'


def getObs(opTable,cID,obs,asNumpy=True):
    """ 
    Grab observable obs for configuration cID from opTable. 
    """
    resl=opTable[cID]['l'][obs]
    ress=opTable[cID]['s'][obs]
    if len(resl) != len(ress): 
        logger.warn(f"len({obs}l) != len({obs}s) cID = "+cID+"... skipping")
        raise ToolboxException
    if asNumpy:
        return np.array(resl), np.array(ress)
    return resl, ress


def _parseOperator(mass,lp,ReOP,ImOP,lVec,sVec,lineno,densFile):
    if rel_check(mass, lp.ml):
        lVec.append(lp.Nc*complex(ReOP, ImOP))
    elif rel_check(mass, lp.ms):
        sVec.append(lp.Nc*complex(ReOP, ImOP))
    else:
        logger.TBRaise("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)


def loadDens(densFile,confID,lp,inTable=None) -> dict:
    """ 
    Allows reading of output from C. Schmidt's Dense code. The Dense code produces measurements of various operators
    relevant for calculating conserved charge fluctuations. We store as a dictionary indexed by confID, which lets us
    conveniently combine data when there are multiple dense files per configuration. Here we update the table
    inTable for the new read-in, which yields outTable. Only supports Nf=2+1 configurations for the time being.

    Parameters
    ----------
    densFile : str
        Name of the Dense code file to be read.
    confID : str
        A unique identifier for the configuration on which dense calculation was performed. Every person seems to have a
        different scheme for labelling configurations, so this needs to be a string to be as flexible as possible.
    lp : latticeParams
        Parameters for the ensemle the configuration belongs to.
    inTable : dict
        A table indexed by confID. Its values are a list of operators that have been measured.

    Returns
    -------
    outTable : dict
        A table indexed by confID. Its values are a list of operators that have been measured.
    """

    if len(confID) != len(confID.strip()):
        logger.TBRaise('confID must not contain whitespace.')

    if inTable is None:
        outTable  = {}
    else:
        checkType(dict,inTable=inTable)
        outTable = inTable

    try:
        infile = open(densFile, 'r')
    except Exception as e:
        logger.warn("Unable to open file",densFile)
        raise e

    logger.info("Loading file",densFile)

    # In the following light and strange quarks are indexed by l and s, respectively.
    # Quark flavor label is suppressed. dM = dM/dmu. We express these operators both
    # in terms of the massive Dirac matrix M and the operators A intoduced in the
    # electric-charge cumulants paper.
    if not confID in outTable:
        trMlVec        = []  # 1  : tr M^-1
        trMsVec        = []
        trMdMlVec      = []  # 2  : tr M^-1 dM                 : tr A
        trMdMsVec      = []
        trMd2MlVec     = []  # 3  : tr M^-1 ddM                : tr A2
        trMd2MsVec     = []
        trMd3MlVec     = []  # 4  : tr M^-1 dddM               : tr A3
        trMd3MsVec     = []
        trMd4MlVec     = []  # 5  : tr M^-1 ddddM              : tr A4
        trMd4MsVec     = []
        trM2lVec       = []  # 6  : tr M^-1 M^-1
        trM2sVec       = []
        trMdMMlVec     = []  # 11 : tr M^-1 d M M^-1  
        trMdMMsVec     = []
        trMdM2lVec     = []  # 12 : tr ( M^-1 dM )^2           : tr A^2
        trMdM2sVec     = []
        trMdMMd2MlVec  = []  # 13 : tr M^-1 dM M^-1 ddM        : tr A A2
        trMdMMd2MsVec  = []
        trMdMMd3MlVec  = []  # 14 : tr M^-1 dM M^-1 dddM       : tr A A3
        trMdMMd3MsVec  = []
        trMdMMd4MlVec  = []  # 15 : tr M^-1 dM M^-1 ddddM      : tr A A4
        trMdMMd4MsVec  = []
        trMd2MMlVec    = []  # 16 : tr M^-1 ddM M^-1 
        trMd2MMsVec    = []
        trMd2MMdMlVec  = []  # 17 : tr M^-1 ddM M^-1 dM        : tr A2 A
        trMd2MMdMsVec  = []
        trMd2M2lVec    = []  # 18 : tr ( M^-1 ddM )^2          : tr A2^2
        trMd2M2sVec    = []
        trMdM2MlVec    = []  # 61 : tr ( M^-1 dM )^2 M^-1 
        trMdM2MsVec    = []
        trMdM3lVec     = []  # 62 : tr ( M^-1 dM )^3           : tr A^3
        trMdM3sVec     = []
        trMdM3Md2MlVec = []  # 63 : tr ( M^-1 dM )^2 M^-1 ddM  : tr A^2 A2
        trMdM3Md2MsVec = []
        trMdM3MlVec    = []  # 311: tr ( M^-1 dM )^3 M^-1 
        trMdM3MsVec    = []
        trMdM4lVec     = []  # 312: tr ( M^-1 dM )^4           : tr A^4
        trMdM4sVec     = []
        
    else:
        trMlVec       , trMsVec        = getObs(outTable,confID,'trM'       ,False) 
        trMdMlVec     , trMdMsVec      = getObs(outTable,confID,'trMdM'     ,False) 
        trMd2MlVec    , trMd2MsVec     = getObs(outTable,confID,'trMd2M'    ,False) 
        trMd3MlVec    , trMd3MsVec     = getObs(outTable,confID,'trMd3M'    ,False)
        trMd4MlVec    , trMd4MsVec     = getObs(outTable,confID,'trMd4M'    ,False)
        trM2lVec      , trM2sVec       = getObs(outTable,confID,'trM2'      ,False)
        trMdMMlVec    , trMdMMsVec     = getObs(outTable,confID,'trMdMM'    ,False)
        trMdM2lVec    , trMdM2sVec     = getObs(outTable,confID,'trMdM2'    ,False) 
        trMdMMd2MlVec , trMdMMd2MsVec  = getObs(outTable,confID,'trMdMMd2M' ,False)
        trMdMMd3MlVec , trMdMMd3MsVec  = getObs(outTable,confID,'trMdMMd3M' ,False)
        trMdMMd4MlVec , trMdMMd4MsVec  = getObs(outTable,confID,'trMdMMd4M' ,False)
        trMd2MMlVec   , trMd2MMsVec    = getObs(outTable,confID,'trMd2MM'   ,False)
        trMd2MMdMlVec , trMd2MMdMsVec  = getObs(outTable,confID,'trMd2MMdM' ,False)
        trMd2M2lVec   , trMd2M2sVec    = getObs(outTable,confID,'trMd2M2'   ,False)
        trMdM2MlVec   , trMdM2MsVec    = getObs(outTable,confID,'trMdM2M'   ,False)
        trMdM3lVec    , trMdM3sVec     = getObs(outTable,confID,'trMdM3'    ,False)
        trMdM3Md2MlVec, trMdM3Md2MsVec = getObs(outTable,confID,'trMdM3Md2M',False)
        trMdM3MlVec   , trMdM3MsVec    = getObs(outTable,confID,'trMdM3M'   ,False)
        trMdM4lVec    , trMdM4sVec     = getObs(outTable,confID,'trMdM4'    ,False)

    # Read in data from the dense file.
    lineno = 0
    for line in infile:

        if 'COMPLETE' in line:
            continue

        lineno += 1
        col = line.split()

        # Parse the densFile.
        try:
            OPID = int(  col[0])
            mass = float(col[2])
            ReOP = float(col[3])
            ImOP = float(col[4])
        except Exception as e:
            infile.close()
            logger.warn("Read error on line", lineno, "of file", densFile)
            raise e

        # If you have a trace squared, you must normalize by volume. Each trace must be 
        # scaled by Nc. Naming convention:
        #     trM  = tr M^-1
        #   trMdM  = tr M^-1 dM/dmu
        #   trMd2M = tr M^-1 d^2M/dmu^2
        # l: light, s: strange
        if OPID == 1:
            _parseOperator(mass,lp,ReOP,ImOP,trMlVec       ,trMsVec       ,lineno,densFile)
        elif OPID == 2:
            _parseOperator(mass,lp,ReOP,ImOP,trMdMlVec     ,trMdMsVec     ,lineno,densFile)
        elif OPID == 3:
            _parseOperator(mass,lp,ReOP,ImOP,trMd2MlVec    ,trMd2MsVec    ,lineno,densFile)
        elif OPID == 4:
            _parseOperator(mass,lp,ReOP,ImOP,trMd3MlVec    ,trMd3MsVec    ,lineno,densFile)  
        elif OPID == 5:
            _parseOperator(mass,lp,ReOP,ImOP,trMd4MlVec    ,trMd4MsVec    ,lineno,densFile) 
        elif OPID == 6:
            _parseOperator(mass,lp,ReOP,ImOP,trM2lVec      ,trM2sVec      ,lineno,densFile) 
        elif OPID == 11:
            _parseOperator(mass,lp,ReOP,ImOP,trMdMMlVec    ,trMdMMsVec    ,lineno,densFile) 
        elif OPID == 12:
            _parseOperator(mass,lp,ReOP,ImOP,trMdM2lVec    ,trMdM2sVec    ,lineno,densFile) 
        elif OPID == 13:
            _parseOperator(mass,lp,ReOP,ImOP,trMdMMd2MlVec ,trMdMMd2MsVec ,lineno,densFile) 
        elif OPID == 14:
            _parseOperator(mass,lp,ReOP,ImOP,trMdMMd3MlVec ,trMdMMd3MsVec ,lineno,densFile) 
        elif OPID == 16:
            _parseOperator(mass,lp,ReOP,ImOP,trMdMMd4MlVec ,trMdMMd4MsVec ,lineno,densFile) 
        elif OPID == 17:
            _parseOperator(mass,lp,ReOP,ImOP,trMd2MMlVec   ,trMd2MMsVec   ,lineno,densFile) 
        elif OPID == 18:
            _parseOperator(mass,lp,ReOP,ImOP,trMd2MMdMlVec ,trMd2MMdMsVec ,lineno,densFile) 
        elif OPID == 61:
            _parseOperator(mass,lp,ReOP,ImOP,trMd2M2lVec   ,trMd2M2sVec   ,lineno,densFile) 
        elif OPID == 62:
            _parseOperator(mass,lp,ReOP,ImOP,trMdM2MlVec   ,trMdM2MsVec   ,lineno,densFile) 
        elif OPID == 63:
            _parseOperator(mass,lp,ReOP,ImOP,trMdM3lVec    ,trMdM3sVec    ,lineno,densFile) 
        elif OPID == 311:
            _parseOperator(mass,lp,ReOP,ImOP,trMdM3Md2MlVec,trMdM3Md2MsVec,lineno,densFile) 
        elif OPID == 312:
            _parseOperator(mass,lp,ReOP,ImOP,trMdM3MlVec   ,trMdM3MsVec   ,lineno,densFile) 
        else:
            continue

    outTable[confID] = { 
        'l': {
            'trM'       : trMlVec, 
            'trMdM'     : trMdMlVec, 
            'trMd2M'    : trMd2MlVec,
            'trMd3M'    : trMd3MlVec,
            'trMd4M'    : trMd4MlVec,
            'trM2'      : trM2lVec,
            'trMdMM'    : trMdMMlVec,
            'trMdM2'    : trMdM2lVec,
            'trMdMMd2M' : trMdMMd2MlVec,
            'trMdMMd3M' : trMdMMd3MlVec,
            'trMdMMd4M' : trMdMMd4MlVec,
            'trMd2MM'   : trMd2MMlVec, 
            'trMd2MMdM' : trMd2MMdMlVec,
            'trMd2M2'   : trMd2M2lVec,  
            'trMdM2M'   : trMdM2MlVec,  
            'trMdM3'    : trMdM3lVec,   
            'trMdM3Md2M': trMdM3Md2MlVec,
            'trMdM3M'   : trMdM3MlVec,
            'trMdM4'    : trMdM4lVec,
            },
        's': { 
            'trM'       : trMsVec, 
            'trMdM'     : trMdMsVec, 
            'trMd2M'    : trMd2MsVec,
            'trMd2M'    : trMd2MsVec,
            'trMd3M'    : trMd3MsVec,
            'trMd4M'    : trMd4MsVec,
            'trM2'      : trM2sVec,
            'trMdMM'    : trMdMMsVec,
            'trMdM2'    : trMdM2sVec,
            'trMdMMd2M' : trMdMMd2MsVec,
            'trMdMMd3M' : trMdMMd3MsVec,
            'trMdMMd4M' : trMdMMd4MsVec,
            'trMd2MM'   : trMd2MMsVec, 
            'trMd2MMdM' : trMd2MMdMsVec,
            'trMd2M2'   : trMd2M2sVec,  
            'trMdM2M'   : trMdM2MsVec,  
            'trMdM3'    : trMdM3sVec,   
            'trMdM3Md2M': trMdM3Md2MsVec,
            'trMdM3M'   : trMdM3MsVec,
            'trMdM4'    : trMdM4sVec,
        }
    }

    infile.close()

    return outTable
