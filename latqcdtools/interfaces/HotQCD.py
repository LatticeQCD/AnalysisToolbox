# 
# HotQCD.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the HotQCD collaboration.
#

import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.math.math import rel_check
from latqcdtools.interfaces.interfaces import HotQCD_MILC_Params


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
    Implementing this as a function makes sure everyone using the Toolbox as the same convention and, more importantly,
    ensures that the tags have no whitespace in them, which otherwise can throw off the column counting of methods in
    the denseObs module. 
    """
    return f'{stream}:{conf}'


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
    Nc = lp.Nc

    # In the following light and strange quarks are indexed by l and s, respectively.
    if not confID in outTable:
        nlVec      = []  # tr M^-1 d M 
        nsVec      = []
        nl2Vec     = []  # tr ( M^-1 d M )^2
        ns2Vec     = []
        MddMlVec   = []  # tr M^-1 dd M
        MddMsVec   = []
        trMinvlVec = []  # tr M^-1
        trMinvsVec = []
    else:
        nlVec      = outTable[confID][0]
        nsVec      = outTable[confID][1]
        nl2Vec     = outTable[confID][2]
        ns2Vec     = outTable[confID][3]
        MddMlVec   = outTable[confID][4]
        MddMsVec   = outTable[confID][5]
        trMinvlVec = outTable[confID][6]
        trMinvsVec = outTable[confID][7]

    # Read in data from the dense file.
    lineno = 0
    for line in infile:

        lineno += 1
        col = line.split()

        # Parse the densFile.
        try:
            OPID = int(col[0])
            mass = float(col[2])
            ReOP = float(col[3])
            ImOP = float(col[4])
        except Exception as e:
            infile.close()
            logger.warn("Read error on line", lineno, "of file", densFile)
            raise e

        # If you have a trace squared, you must normalize by volume. Each trace must be normalized by Nc.
        if OPID == 1:  # tr M^-1
            if rel_check(mass, lp.ml):
                trMinvlVec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                trMinvsVec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBRaise("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        elif OPID == 2:  # D1 = tr M^-1 d M
            if rel_check(mass, lp.ml):
                nlVec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                nsVec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBRaise("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        elif OPID == 3:  # tr M^-1 dd M
            if rel_check(mass, lp.ml):
                MddMlVec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                MddMsVec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBRaise("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        elif OPID == 12:  # tr (M^-1 d M)**2
            if rel_check(mass, lp.ml):
                nl2Vec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                ns2Vec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBRaise("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        else:
            continue

    outTable[confID] = [nlVec, nsVec, nl2Vec, ns2Vec, MddMlVec, MddMsVec, trMinvlVec, trMinvsVec]

    infile.close()

    return outTable
