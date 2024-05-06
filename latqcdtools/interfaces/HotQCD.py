# 
# HotQCD.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the HotQCD collaboration.
#
import latqcdtools.base.logger as logger
from latqcdtools.math.math import rel_check
from latqcdtools.interfaces.interfaces import HotQCD_MILC_Params


class HotQCDParams(HotQCD_MILC_Params):
    """A class to handle and check the input parameters of a lattice run, especially for HotQCD."""

    def __repr__(self) -> str:
        return "HotQCDParams"

    def getcGradFlowPureGauge(self) -> str:
        return 's'+str(self.Ns).zfill(3)+'t'+str(self.Nt).zfill(2)+'_b0'+self.cbeta+'00'


def massRatioToMasses(msml, Nt, cbeta, Nf='21'):
    """ This is a way to get either the light+strange masses or quark+preconditioner masses
        given the ratio ms/ml or mpre/mf. This uses HotQCD ensembles, listed in tables below.

    Args:
        msml (int)
        Nt (int)
        cbeta (str): Beta as a character. e.g. 6.345 is 6345 
        Nf (str): Number of dynamical fermions. Defaults to '21'.

    Returns:
        str, str: either (ml, ms) or (mf, mpre) 
    """
    if msml is None:
        cm1=None
        cm2=None
    else:
        massTable=quarkMassTableHISQ(Nf,Nt,msml)
        try:    
            cm1=massTable[cbeta][0]
            cm2=massTable[cbeta][1]
        except KeyError:
            _badbeta(cbeta,msml,Nt,Nf)
    return cm1, cm2


def _badbeta(beta,msml,Nt,Nf):
    logger.TBError("No entries for beta =",beta,", Nf =",Nf,", Nt =",Nt,", msml =",msml)
def _badmsml(msml,Nt,Nf):
    logger.TBError("No entries for Nf =",Nf,", Nt =",Nt,", msml =",msml)
def _badNt(Nt,Nf):
    logger.TBError("No entries for Nf =",Nf,", Nt =",Nt)
def _badNf(Nf):
    logger.TBError("No entries for Nf =",Nf)


def quarkMassTableHISQ(Nf, Nt, msml) -> dict:
    """Lookup tables for HotQCD HISQ runs."""
    Table = None

    if Nf=='21':

        if Nt==12:

            if msml==27:
                Table = {'6794': ['00167', '0450'],
                         '6850': ['00157', '0424'],
                         '6910': ['00148', '0401']}
            else:
                _badmsml(msml,Nt,Nf)

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
                         '6500': ['00228', '0614'],
                         '6640': ['00196', '0528']}
            else:
                _badmsml(msml,Nt,Nf)

        elif Nt==6:
            if msml==27:
                Table = {'6325': ['00277' ,'0748'],
                         '6245': ['00307' ,'0830'],
                         '5980': ['00435' ,'1173'],
                         '6300': ['00285' ,'0770'],
                         '6170': ['00336' ,'0908'],
                         '6120': ['00357' ,'0963'], # a guess based on a spline interpolation
                         '6038': ['00399' ,'1077'], # a guess based on a spline interpolation
                         '6274': ['00295' ,'0798'], # a guess based on a spline interpolation
                         '5850': ['005274','1424'],
                }
            else:
                _badmsml(msml,Nt,Nf)
        
        elif Nt==4:
            if msml==27:
                Table = {'5850': ['005274', '1424'], 
                         '5900': ['004889', '1320'],
                         '5925': ['004722', '1275'],
                         '5950': ['004556', '1230'],
                         '5960': ['004486', '1211'],
                         '5970': ['004416', '1192'],
                         '5973': ['004394', '1186'], # a guess based on a spline interpolation
                         '5975': ['004377', '1183'],
                         '5980': ['004346', '1173'],
                         '5990': ['004279', '1155'],
                         '6000': ['004215', '1138'],
                         '6025': ['004074', '1100']} 
            else:
                _badmsml(msml,Nt,Nf)

        else:
            _badNt(Nt,Nf) 

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
                _badmsml(msml,Nt,Nf)

        elif Nt==16:

            if msml==27:
                Table = {'6050': ['00394', '1064'],
                         '6315': ['00281', '07587']}
            else:
                _badmsml(msml,Nt,Nf)

        else:
            _badNt(Nt,Nf) 

    elif Nf=='5':

        if Nt==6:

            if msml==25:
                Table = {'4525': ['002', '05'],
                         '4550': ['002', '05'],
                         '4575': ['002', '05'],
                         '4600': ['002', '05'],
                         '4610': ['002', '05'],
                         '4620': ['002', '05'],
                         '4625': ['002', '05'],
                         '4630': ['002', '05'],
                         '4632': ['002', '05'],
                         '4635': ['002', '05'],
                         '4637': ['002', '05'],
                         '4640': ['002', '05'],
                         '4650': ['002', '05'],
                         '4660': ['002', '05'],
                         '4670': ['002', '05'],
                         '4675': ['002', '05'],
                         '4680': ['002', '05'],
                         '4690': ['002', '05'],
                         '4700': ['002', '05']} 
            else:
                _badmsml(msml,Nt,Nf)

        else:
            _badNt(Nt,Nf) 

    else:
        _badNf(Nf) 

    return Table


def makeConfTag(conf,stream) -> str:
    """ This takes a configuration number conf and stream label stream to make a tag labelling a configuration.
    Implementing this as a function makes sure everyone using the Toolbox as the same convention and, more importantly,
    ensures that the tags have no whitespace in them, which otherwise can throw off the column counting of methods in
    the denseObs module. """
    return str(stream)+':'+str(conf)


def loadDens(densFile,confID,lp,inTable=None) -> dict:
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

    if len(confID) != len(confID.strip()):
        logger.TBError('confID must not contain whitespace.')

    logger.warn('This may be wrong. Do not use for now.')

    if inTable is None:
        outTable  = {}
    elif not isinstance(inTable,dict):
        logger.TBError("Must pass dict to inTable, or else pass None.")
    else:
        outTable = inTable

    try:
        infile = open(densFile, 'r')
    except Exception as e:
        logger.warn("Unable to open file",densFile)
        raise e

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
