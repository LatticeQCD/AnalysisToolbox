# 
# testDensObs.py
# 
# D. Clarke
# 
# Test to see if we can read in some dense files and compute observables from that.
# 

from latqcdtools.interfaces.HotQCD import HotQCDParams, loadDens
from latqcdtools.physics.denseObs import op_to_obs
import latqcdtools.base.logger as logger
from latqcdtools.base.cleanData import restrictAtCol
from latqcdtools.math.math import rel_check
from latqcdtools.base.readWrite import readTable
from latqcdtools.testing import concludeTest
from latqcdtools.base.utilities import unvector, ls

logger.set_log_level('INFO')

lp = HotQCDParams(Nsigma=40, Ntau=8, coupling='6260',mass1='002025',mass2='0810')

EPSILON = 1e-5

def testDensObs():

    initialize = True

    for stream in [0, 1]:

        directory = f'denseObs/str{stream}'

        for filename in ls(f'{directory}/*'):

            confID = filename.split('_')[1]

            if initialize:
                opTable = loadDens(filename, confID, lp)
                initialize = False
            else:
                opTable = loadDens(filename, confID, lp, opTable)

    OBS = op_to_obs(opTable, lp, writeFiles=False)

    REFdata = readTable('denseObs/n_n2_dn_table_ms40_b6260.d')


    lpass = True
    for iconf in range(len(OBS['confID'])):

        confID = OBS['confID'][iconf]

        stream = int(confID.split('.')[0])
        conf   = int(confID.split('.')[1])

        REFdata_sc = restrictAtCol( restrictAtCol(REFdata,0,stream), 1, conf )

        REFchi2l   = unvector( REFdata_sc[8]   )*lp.Nt**2
        REFchi2s   = unvector( REFdata_sc[9]   )*lp.Nt**2 
        REFchi11ll = unvector( REFdata_sc[10]  )*lp.Nt**2 
        REFchi11ls = unvector( REFdata_sc[11]  )*lp.Nt**2 
        REFchi2B   = unvector( REFdata_sc[12]  )*lp.Nt**2 
        REFchi2Q   = unvector( REFdata_sc[13]  )*lp.Nt**2 

        chi2l   = OBS['X2l'][iconf]
        chi2s   = OBS['X2s'][iconf]
        chi11ll = OBS['X11ll'][iconf]
        chi11ls = OBS['X11ls'][iconf]
        chi2B   = OBS['X2B'][iconf]
        chi2Q   = OBS['X2Q'][iconf]

        if not rel_check(chi2l,REFchi2l):
            lpass = False
            logger.TBFail('chi2l',chi2l,'ref',REFchi2l,'conf',confID)

        if not rel_check(chi2s,REFchi2s):
            lpass = False
            logger.TBFail('chi2s',chi2s,'ref',REFchi2s,'conf',confID)

        if not rel_check(chi11ll,REFchi11ll):
            lpass = False
            logger.TBFail('chi11ll',chi11ll,'ref',REFchi11ll,'conf',confID)

        if not rel_check(chi11ls,REFchi11ls):
            lpass = False
            logger.TBFail('chi11ls',chi11ls,'ref',REFchi11ls,'conf',confID)

        if not rel_check(chi2B,REFchi2B):
            lpass = False
            logger.TBFail('chi2B',chi2B,'ref',REFchi2B,'conf',confID)

        if not rel_check(chi2Q,REFchi2Q):
            lpass = False
            logger.TBFail('chi2Q',chi2Q,'ref',REFchi2Q,'conf',confID)

    concludeTest(lpass)


if __name__ == '__main__':
    testDensObs()