# 
# HotQCD.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the HotQCD collaboration.
#
import latqcdtools.base.logger as logger
from latqcdtools.physics.lattice_params import latticeParams


class HotQCDParams(latticeParams):
    """A class to handle and check the input parameters of a lattice run, especially for HotQCD."""

    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2
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

    if Nf=='21':

        if Nt==12:

            if msml==27:
                Table = {'6794': ['00167', '0450'],
                         '6850': ['00157', '0424'],
                         '6910': ['00148', '0401']}
            else:
                logger.TBError("Invalid ms/ml for quark mass table.")

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
                logger.TBError("Invalid ms/ml for quark mass table.")

        else:
            logger.TBError("Invalid Nt for quark mass table.")

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
                logger.TBError("Invalid ms/ml for quark mass table.")

        elif Nt==16:

            if msml==27:
                Table = {'6050': ['00394', '1064'],
                         '6315': ['00281', '07587']}
            else:
                logger.TBError("Invalid ms/ml for quark mass table.")

        else:
            logger.TBError("Invalid Nt for quark mass table.")

    elif Nf=='5':

        if Nt==6:

            if msml==25:
                Table = {'4637': ['002', '05'] }
            else:
                logger.TBError("Invalid ms/ml for quark mass table.")

        else:
            logger.TBError("Invalid Nt for quark mass table.")


    else:
        logger.TBError("Invalid Nf for quark mass table.")

    return Table