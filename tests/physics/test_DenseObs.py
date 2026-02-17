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
from latqcdtools.base.utilities import unvector, ls, toNumpy
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.base.speedify import parallel_function_eval
import numpy as np

logger.set_log_level('INFO')

lp = HotQCDParams(Nsigma=40, Ntau=8, coupling='6260',mass1='002025',mass2='0810')

EPSILON = 1e-5

def compare(test,REF,name,cID) -> bool:
    if not rel_check(test,REF):
        logger.TBFail(f'{name} {test} ref {REF} conf {cID}')
        return False
    return True

def compareJack(arr,REF,REFe,name) -> bool:
    passed=True
    m, e = jackknife(np.mean,arr,numb_blocks=20)
    if not rel_check(m.real,REF,prec=3e-2):
        passed = False
        logger.TBFail(f'{name} {m.real} ref {REF}')
    if not rel_check(m.real,REF,prec=3e-2):
        passed = False
        logger.TBFail(f'{name}e {e} refe {REFe}')
    return passed

def constructOpTable(opTableList,streamList) -> dict:
    opTable = opTableList[0]
    for stream in streamList[1:]:
        for conf in opTableList[stream]:
            if not conf in opTable:
                opTable[conf] = {}
            for q in ['l','s']:
                if not q in opTable[conf]:
                    opTable[conf][q]= {}
                for obs in opTableList[stream][conf][q]:
                    opTable[conf][q][obs]=opTableList[stream][conf][q][obs]
    return opTable


def testDensObs():

    def processDense(stream):

        directory = f'denseObs/str{stream}'

        initialize = True
        for filename in ls(f'{directory}/*'):
            confID = filename.split('_')[1]
            if initialize:
                opTable = loadDens(filename, confID, lp)
                initialize = False
            else:
                opTable = loadDens(filename, confID, lp, opTable)

        return opTable

    temp = parallel_function_eval(processDense,[0,1]) 
    opTable = constructOpTable(temp,[0,1]) 

    OBS = op_to_obs(opTable, lp, writeFiles=False)

    REFdata = readTable('denseObs/n_n2_dn_table_ms40_b6260.d')


    lpass = True

    chi2lavg,chi2savg,chi11llavg,chi11lsavg,chi2Bavg,chi2Qavg = [],[],[],[],[],[]

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

        chi2lavg.append(chi2l)  
        chi2savg.append(chi2s)  
        chi11llavg.append(chi11ll) 
        chi11lsavg.append(chi11ls) 
        chi2Bavg.append(chi2B) 
        chi2Qavg.append(chi2Q) 

        lpass *= compare(chi2l  ,REFchi2l  ,'chi2l'  ,confID)
        lpass *= compare(chi2s  ,REFchi2s  ,'chi2s'  ,confID)
        lpass *= compare(chi11ll,REFchi11ll,'chi11ll',confID)
        lpass *= compare(chi11ls,REFchi11ls,'chi11ls',confID)
        lpass *= compare(chi2B  ,REFchi2B  ,'chi2B'  ,confID)
        lpass *= compare(chi2Q  ,REFchi2Q  ,'chi2Q'  ,confID)

    chi2lavg,chi2savg,chi11llavg,chi11lsavg,chi2Bavg,chi2Qavg = toNumpy( chi2lavg,chi2savg,chi11llavg,chi11lsavg,chi2Bavg,chi2Qavg ) 

    # 
    # Comparison with Jishnu's code
    #
    lpass *= compareJack(chi2lavg  , 3.02582962e-01,2.39426453e-02,'chi2l'  )
    lpass *= compareJack(chi2savg  , 1.10486541e-01,3.64726294e-03,'chi2s'  )
    lpass *= compareJack(chi11llavg,-6.20992517e-02,2.24982049e-02,'chi11ll')
    lpass *= compareJack(chi2Bavg  , 5.75703687e-02,1.30232279e-02,'chi2B'  )
    lpass *= compareJack(chi2Qavg  , 2.12050964e-01,7.11385988e-03,'chi2Q'  )

    concludeTest(lpass)


if __name__ == '__main__':
    testDensObs()