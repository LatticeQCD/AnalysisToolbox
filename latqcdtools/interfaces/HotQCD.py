# 
# HotQCD.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the HotQCD collaboration.
#
import latqcdtools.base.logger as logger
from latqcdtools.base.check import rel_check
from latqcdtools.interfaces.interfaces import HotQCD_MILC_Params


class HotQCDParams(HotQCD_MILC_Params):
    """A class to handle and check the input parameters of a lattice run, especially for HotQCD."""

    def getcGradFlowPureGauge(self):
        return 's'+str(self.Ns).zfill(3)+'t'+str(self.Nt).zfill(2)+'_b0'+self.cbeta+'00'


def massRatioToMasses(msml, Nt, cbeta, Nf='21'):
    """Get mass parameters given the others."""

    if msml is None:
        cm1=None
        cm2=None
    else:
        massTable=quarkMassTableHISQ(Nf,Nt,msml)
        cm1=massTable[cbeta][0]
        cm2=massTable[cbeta][1]
    return cm1, cm2


def quarkMassTableHISQ(Nf, Nt, msml):
    """Lookup tables for HotQCD HISQ runs."""
    Table = None

    if Nf=='21':

        if Nt==12:

            if msml==27:
                Table = {'6794': ['00167', '0450'],
                         '6850': ['00157', '0424'],
                         '6910': ['00148', '0401']}
            else:
                logger.TBError("Invalid ms/ml.")

        elif Nt==8:

            if msml==160:
                Table = {'6285': ['00049375', '0790'],
                         '6300': ['00048250', '0772'],
                         '6315': ['00047500', '0760'],
                         '6330': ['00046625', '0746'],
                         '6354': ['00045500', '0728'],
                         '6372': ['00044456', '0711'],
                         '6390': ['00043375', '0694'],
                         '6423': ['00041875', '0670'],
                         '6445': ['00040750', '0652']}
            elif msml==80:
                Table = {'6285': ['0009875', '0790'],
                         '6300': ['0009650', '0772'],
                         '6315': ['0009500', '0760'],
                         '6330': ['0009325', '0746'],
                         '6354': ['0009100', '0728'],
                         '6372': ['0008891', '0711'],
                         '6390': ['0008675', '0694'],
                         '6423': ['0008375', '0670'],
                         '6445': ['0008150', '0652']}
            elif msml==40:
                Table = {'6260': ['002025', '0810'],
                         '6285': ['001975', '0790'],
                         '6300': ['001930', '0772'],
                         '6315': ['001900', '0760'],
                         '6330': ['001865', '0746'],
                         '6354': ['001820', '0728'],
                         '6365': ['001790', '0716'],
                         '6390': ['001735', '0694'],
                         '6423': ['001675', '0670'],
                         '6445': ['001630', '0652'],
                         '6474': ['001580', '0632'],
                         '6500': ['001535', '0614']}
            elif msml==27:
                Table = {'6245': ['00307', '0830'],
                         '6285': ['00293', '0790'],
                         '6315': ['00281', '0759'],
                         '6354': ['00270', '0728'],
                         '6390': ['00257', '0694'],
                         '6423': ['00248', '0670'],
                         '6445': ['00241', '0652'],
                         '6474': ['00234', '0632'],
                         '6500': ['00228', '0614']}
            else:
                logger.TBError("Invalid ms/ml.")

        else:
            logger.TBError("Invalid Nt.")

    elif Nf=='3':

        if Nt==8:

            if msml==27:
                Table = {'6050': ['00394', '1064'],
                         '6085': ['00376', '1015'],
                         '6110': ['00364', '0982'],
                         '6160': ['00341', '0920'],
                         '6175': ['00334', '0902'],
                         '6190': ['00328', '0885'],
                         '6215': ['00318', '0858'],
                         '6245': ['00307', '08289'],
                         '6285': ['00293', '07911'],
                         '6315': ['00281', '07587']}
            else:
                logger.TBError("Invalid ms/ml.")

        elif Nt==16:

            if msml==27:
                Table = {'6050': ['00394', '1064'],
                         '6315': ['00281', '07587']}
            else:
                logger.TBError("Invalid ms/ml.")

        else:
            logger.TBError("Invalid Nt.")

    elif Nf=='5':

        if Nt==6:

            if msml==25:
                Table = {'4637': ['002', '05'] }
            else:
                logger.TBError("Invalid ms/ml.")

        else:
            logger.TBError("Invalid Nt.")


    else:
        logger.TBError("Invalid Nf.")

    return Table


def loadDens(densFile,confID,lp,inTable=None):
    """ Allows reading of output from C. Schmidt's Dense code. The Dense code produces measurements of various operators
    relevant calculating conserved charge fluctuations. We store as a dictionary indexed by confID, which lets us
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

    if inTable is None:
        outTable  = {}
    elif not isinstance(inTable,dict):
        logger.TBError("Must pass dict to inTable, or else pass None.")
    else:
        outTable = inTable

    try:
        infile = open(densFile, 'r')
    except Exception as e:
        logger.TBError("Unable to open file",densFile)

    Nc = lp.Nc

    # In the following light and strange quarks are indexed by l and s, respectively.
    if not confID in outTable:
        nlVec      = []  # quark number density, n
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
        except ValueError:
            logger.TBError("Read error on line", lineno, "of file", densFile)
        try:
            mass = float(col[2])
            ReOP = float(col[3])
            ImOP = float(col[4])
        except IndexError:
            logger.TBError("Read error on line", lineno, "of file", densFile)

        # If you have a trace squared, you must normalize by volume. Each trace must be normalized by Nc.
        if OPID == 1:  # tr M^-1
            if rel_check(mass, lp.ml):
                trMinvlVec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                trMinvsVec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBError("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        elif OPID == 2:  # D1 = tr M^-1 d M
            if rel_check(mass, lp.ml):
                nlVec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                nsVec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBError("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        elif OPID == 3:  # tr M^-1 dd M
            if rel_check(mass, lp.ml):
                MddMlVec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                MddMsVec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBError("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        elif OPID == 12:  # tr (M^-1 d M)**2
            if rel_check(mass, lp.ml):
                nl2Vec.append(Nc * complex(ReOP, ImOP))
            elif rel_check(mass, lp.ms):
                ns2Vec.append(Nc * complex(ReOP, ImOP))
            else:
                logger.TBError("Unexpected mass on line", lineno, "of file", densFile, ". ms, ml, m =", lp.ms, lp.ml, mass)

        else:
            continue

        outTable[confID] = [nlVec, nsVec, nl2Vec, ns2Vec, MddMlVec, MddMsVec, trMinvlVec, trMinvsVec]

    infile.close()

    return outTable
