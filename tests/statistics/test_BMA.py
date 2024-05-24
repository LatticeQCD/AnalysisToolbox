# 
# testBMA.py                                                               
# 
# D. Clarke 
# 
# Quick test of BMA. Test by applying to practical example IIIC of
# 10.1103/PhysRevD.103.114502. Implementation from the authors can
# be found here: https://github.com/etneil/model_average_paper.
# 


import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.printErrorBars import getValuesFromErrStr
from latqcdtools.base.initialize import DEFAULTSEED
from latqcdtools.base.plotting import plt, plot_dots, clearPlot
from latqcdtools.statistics.statistics import plot_func, getModelWeights, modelAverage
from latqcdtools.statistics.fitting import Fitter
from latqcdtools.testing import concludeTest, print_results


logger.set_log_level('INFO')

def F(x):
    """ Model truth """
    return 1.8-0.53*(x/16)+0.31*(x/16)**2

def f_m(x,coeffs):
    """ Model to try. """
    result = 0.
    for j in range(len(coeffs)):
        result += coeffs[j]*(x/16)**j
    return result

paper_data  = ['1.84(13)', '1.90(15)', '1.73(13)', '1.75(12)', '1.53(13)', '1.62(14)', '1.36(13)', 
               '1.67(13)', '1.61(14)', '1.35(13)', '1.71(13)', '1.44(12)', '1.54(12)', '1.62(12)', '1.31(11)']
paper_a0    = ['1.587(32)', '1.803(67)', '1.89(11)', '2.01(16)', '1.98(17)', '1.94(18)']
BAICPAPER   = [30.85, 19.17, 20.23, 20.88, 22.22, 23.79]
PRBAICPAPER = [0.0012, 0.43, 0.25, 0.18, 0.09, 0.04]
BAICCONTROL = [29.21641305840999, 19.229646697947416, 20.126005931878353,
                   20.36733203289231, 20.654635925777384, 22.64892748439832 ]


def testBMA():

    rng = np.random.default_rng(DEFAULTSEED)

    Nt       = 15
    Nmodels  = 6

    SHOWTRUE = False 
    SHOWFITS = False 

    X = np.array(range(1,Nt+1))

    Ym = []
    Ye = []
    for t in range(Nt):
        mean, err = getValuesFromErrStr(paper_data[t])
        Ym.append(mean)
        Ye.append(err)

    if SHOWTRUE:    
        plot_dots(X,Ym,Ye)
        plot_func(F,domain=(0,16))
        plt.show()
        clearPlot()

    fitter = Fitter(f_m,X,Ym,Ye,test_tol=1e-8)

    lpass = True

    PRTEST = getModelWeights(BAICPAPER)
    lpass *= print_results(PRTEST,PRBAICPAPER,prec=1e-1)

    for m in range(Nmodels):
        startingGuess = np.repeat(1,m+1)
        res, res_err, chi_dof, stats=fitter.try_fit(start_params=startingGuess, algorithms=['curve_fit'], detailedInfo=True) 
        BAICTEST = stats['BAIC']
        if (BAICTEST/BAICPAPER[m]) > 1.1:
            lpass = False
        lpass *= print_results(BAICTEST,BAICCONTROL[m],text='BAIC, chi2, m = '+str(m)) 
        if SHOWFITS:
            fitter.plot_data()
            fitter.plot_fit(domain=(0,16))
            plt.show()
            clearPlot()

    am, ae = [], []
    for t in range(Nmodels):
        mean, err = getValuesFromErrStr(paper_a0[t])
        am.append(mean)
        ae.append(err)

    aCONTROL = [1.8847300219753202, 0.14328657401006523]
    aTEST    = modelAverage(am,ae,BAICPAPER)
    lpass *= print_results(aTEST,aCONTROL,text='modelAverage')

    concludeTest(lpass)


if __name__ == '__main__':
    testBMA()