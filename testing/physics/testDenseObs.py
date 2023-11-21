# 
# testDensObs.py
# 
# D. Clarke
# 
# Test to see if we can read in some dense files and compute observables from that.
# 

import glob
import numpy as np
from latqcdtools.interfaces.HotQCD import HotQCDParams, loadDens
from latqcdtools.physics.denseObs import op_to_obs, observablesOfInterest
import latqcdtools.base.logger as logger
from latqcdtools.base.cleanData import restrictAtCol
from latqcdtools.math.math import rel_check
from latqcdtools.base.readWrite import readTable
from latqcdtools.testing import concludeTest
from latqcdtools.base.utilities import unvector
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.statistics.statistics import std_mean, std_err

logger.set_log_level('INFO')

lp = HotQCDParams(Nsigma=40, Ntau=8, coupling='6260',mass1='002025',mass2='0810')

EPSILON = 1e-5

def testDensObs():

    initialize = True

    # Initialize the dense observables you want to measure.
    obs = observablesOfInterest(["confID",
                                 "Nl", "NB", "NQ", "NS", "NI",
                                 "dN/dmul", "dN/dmus", "dN/dmuB", "dN/dmuQ", "dN/dmuS", "dN/dmuI",
                                 "Nl^2", "NB^2", "NQ^2", "NS^2",
                                 "chi2l", "chi2s", "chi11ll", "chi11ls", "chi2B", "chi2Q"])

    logger.info('Testing for mu=0...')
    for stream in [0, 1]:

        directory = 'denseObs/str' + str(stream)

        for filename in glob.iglob(f'{directory}/*'):

            logger.details('Process file', filename)
            confID = filename.split('_')[1]

            if initialize:
                opTable = loadDens(filename, confID, lp)
                initialize = False
            else:
                opTable = loadDens(filename, confID, lp, opTable)

    op_to_obs(opTable, lp, filename='denseObs/denseObservables.d')

    REFdata = readTable('denseObs/n_n2_dn_table_ms40_b6260.d')
    data = np.genfromtxt('denseObs/denseObservables.d', dtype=obs.dtypes, unpack=True)

    # see whether it actually decreases the error to do the imaginary thing.
    Narr = []

    # Compare new class to reference
    lpass = True
    for iconf in range(len(data[0])):

        confID = data[0][iconf]

        stream = int(confID.split('.')[0])
        conf   = int(confID.split('.')[1])

        REFdata_sc = restrictAtCol( restrictAtCol(REFdata,0,stream), 1, conf )

        # When the reference data were calculated, n_B used a different convention, because I was following a paper
        # that did not include a factor 1/3.
        REFReN     = unvector( REFdata_sc[2]/3 )
        REFImN     = unvector( REFdata_sc[3]/3 )
        REFReN2    = unvector( REFdata_sc[4]/9 )
        REFImN2    = unvector( REFdata_sc[5]/9 )
        REFReDN    = unvector( REFdata_sc[6]/9 )
        REFImDN    = unvector( REFdata_sc[7]/9 )
        REFchi2l   = unvector( REFdata_sc[8]   )
        REFchi2s   = unvector( REFdata_sc[9]   )
        REFchi11ll = unvector( REFdata_sc[10]  )
        REFchi11ls = unvector( REFdata_sc[11]  )
        REFchi2B   = unvector( REFdata_sc[12]  )
        REFchi2Q   = unvector( REFdata_sc[13]  )

        ReN     = data[obs.getCol('Re','NB')][iconf]
        ImN     = data[obs.getCol('Im','NB')][iconf]
        ReN2    = data[obs.getCol('Re','NB^2')][iconf]
        ImN2    = data[obs.getCol('Im','NB^2')][iconf]
        ReDN    = data[obs.getCol('Re','dN/dmuB')][iconf]
        ImDN    = data[obs.getCol('Im','dN/dmuB')][iconf]
        chi2l   = data[obs.getCol('Re','chi2l')][iconf]
        chi2s   = data[obs.getCol('Re','chi2s')][iconf]
        chi11ll = data[obs.getCol('Re','chi11ll')][iconf]
        chi11ls = data[obs.getCol('Re','chi11ls')][iconf]
        chi2B   = data[obs.getCol('Re','chi2B')][iconf]
        chi2Q   = data[obs.getCol('Re','chi2Q')][iconf]

        Narr.append(chi2l) 

        if not rel_check(ReN,REFReN,prec=EPSILON):
            lpass = False
            logger.TBFail('Re NB',ReN,'ref',REFReN,'conf',confID)

        if not rel_check(ImN,REFImN,prec=EPSILON):
            lpass = False
            logger.TBFail('Im NB',ImN,'ref',REFImN,'conf',confID)

#        logger.warn('Sign here needs to be checked against imaginary mu data; reference seems wrong.')
#        if not rel_check(ReN2,REFReN2,prec=EPSILON):
#            lpass = False
#            logger.TBFail('Re NB^2',ReN2,'ref',REFReN2,'conf',confID)

        if not rel_check(ImN2,REFImN2,prec=EPSILON):
            lpass = False
            logger.TBFail('Im NB^2',ImN2,'ref',REFImN2,'conf',confID)

        if not rel_check(ReDN,REFReDN,prec=EPSILON):
            lpass = False
            logger.TBFail('Re dN/dmuB',ReDN,'ref',REFReDN,'conf',confID)

        if not rel_check(ImDN,REFImDN,prec=EPSILON):
            lpass = False
            logger.TBFail('Im dN/dmuB',ImDN,'ref',REFImDN,'conf',confID)

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

    print(get_err_str(std_mean(Narr),std_err(Narr))) 

    concludeTest(lpass)


if __name__ == '__main__':
    testDensObs()