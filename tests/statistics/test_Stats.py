# 
# testStats.py                                                               
# 
# D. Clarke
# 
# Tests of some of the basic statistics methods.
#

import numpy as np
import scipy as sp
from latqcdtools.statistics.statistics import gaudif, studif, pearson, std_mean, cov_to_cor, confidence_ellipse,\
    KSTest_1side, KSTest_2side, covariance, std_var, symmetrizeError, forcePositiveSemidefinite
from latqcdtools.math.math import isPositiveSemidefinite
import latqcdtools.base.logger as logger
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.base.plotting import plt
from latqcdtools.base.initialize import DEFAULTSEED
from latqcdtools.base.utilities import toNumpy

logger.set_log_level('INFO')

eps=1e-7

# Some measurements and error
x1=0.4655
e1=0.0088
x2=0.501
e2=0.045
x3=0.480
e3=0.023
x4=0.4745
e4=0.0063

# Some invented numbers of data
ndat1 = 15
ndat2 = 9

# Results produced by software of "Markov Chain Monte Carlo Simulations and Their Statistical Analysis, World
# Scientific, Singapore, 2004.
q12control=0.4387984
q13control=0.5559897
q14control=0.4056413
s12control=0.33726853

# An example data set with CO2 concentrations and global temperature anomaly
ts1 = [284.7,   284.9,   285. ,   285.1,   285.3,   285.4,   285.6,   285.7,   285.9,
       286.1,   286.2,   286.4,   286.5,   286.6,   286.8,   286.9,   287. ,   287.1,
       287.2,   287.4,   287.5,   287.7,   287.9,   288.1,   288.4,   288.7,   289.,
       289.4,   289.8,   290.2,   290.7,   291.2,   291.7,   292.1,   292.6,   293.,
       293.3,   293.6,   293.8,   294. ,   294.2,   294.3,   294.5,   294.6,   294.7,
       294.8,   294.9,   295. ,   295.2,   295.5,   295.8,   296.1,   296.5,   296.8,
       297.2,   297.6,   298.1,   298.5,   298.9,   299.3,   299.7,   300.1,   300.4,
       300.8,   301.1,   301.4,   301.7,   302.1,   302.4,   302.7,   303. ,   303.4,
       303.8,   304.1,   304.5,   305. ,   305.4,   305.8,   306.3,   306.8,   307.2,
       307.7,   308.2,   308.6,   309. ,   309.4,   309.8,   310. ,   310.2,   310.3,
       310.4,   310.4,   310.3,   310.2,   310.1,   310.1,   310.1,   310.2,   310.3,
       310.5,   310.7,   311.1,   311.5,   311.9,   312.4,   313. ,   313.6,   314.2,
       314.9,   315.7,   316.6,   317.3,   318.0,   318.6,   319.4,   320.0,   321.0, 
       321.9,   322.9,   324.2,   325.2,   326.0,   327.1,   328.8,   329.6,   330.7,
       331.7,   333.2,   334.5,   336.8,   338.7,   340.1,   341.4,   343.1,   344.8,
       346.3,   347.6,   349.3,   351.6,   353.2,   354.4,   355.7,   356.5,   357.2,
       358.9,   360.9,   362.7,   363.8,   366.8,   368.5,   369.7,   371.3,   373.4,
       375.9,   377.7,   379.9,   382.0,   384.0,   385.8,   387.6,   390.1,   391.8,
       394.0,   396.7,   398.8,   401.0,   404.4,   406.7,   408.7,   411.6,   414.2,
       416.4,   418.5]
ts2 = [-0.41765878, -0.2333498 ,  -0.22939907,  -0.27035445,  -0.29163003,  -0.2969512,
       -0.32035372, -0.46723005,  -0.3887657 ,  -0.28119546,  -0.39016518,  -0.42927712,
       -0.53639776, -0.3443432 ,  -0.4654367 ,  -0.33258784,  -0.34126064,  -0.35696334,
       -0.35196072, -0.31657043,  -0.32789087,  -0.3685807 ,  -0.32804197,  -0.34133235,
       -0.3732512 , -0.37562594,  -0.42410994,  -0.10110883,  -0.01131519,  -0.30363432,
       -0.31583208, -0.23224552,  -0.29553008,  -0.3464744 ,  -0.49232006,  -0.47112358,
       -0.42090362, -0.49878576,  -0.37937889,  -0.24989556,  -0.50685817,  -0.40131494,
       -0.5075585 , -0.49461925,  -0.48376393,  -0.4487516 ,  -0.28400728,  -0.25980017,
       -0.48579213, -0.35543364,  -0.23447904,  -0.29342857,  -0.43898427,  -0.5333264,
       -0.5975614 , -0.40775132,  -0.3191393 ,  -0.5041577 ,  -0.5138707 ,  -0.5357649,
       -0.5310242 , -0.5392051 ,  -0.47567302,  -0.46715254,  -0.2625924 ,  -0.19184391,
       -0.42020997, -0.54301953,  -0.42458433,  -0.32551822,  -0.2985808 ,  -0.24067703,
       -0.33922812, -0.31793055,  -0.3120622 ,  -0.28242525,  -0.12283547,  -0.22940508,
       -0.20676155, -0.39275664,  -0.1768054 ,  -0.10339768,  -0.14546166,  -0.32234442,
       -0.17433685, -0.20605922,  -0.16952093,  -0.01919893,  -0.01220073,  -0.04079717,
        0.07593584,  0.03812934,   0.00140609,   0.00641407,   0.14410514,   0.04308837,
       -0.1188128 , -0.09120554,  -0.12466127,  -0.14380224,  -0.22662179,  -0.06115397,
        0.01535457,  0.07763074,  -0.11675021,  -0.19730993,  -0.2631656 ,  -0.03533493,
       -0.01763255, -0.04800483,  -0.11548702,  -0.01999739,  -0.06405444,  -0.03680589,
       -0.30586675, -0.2043879 ,  -0.14888458,  -0.11751631,  -0.1686323 ,  -0.03136671,
       -0.08510657, -0.20593274,  -0.0938271 ,   0.04993336,  -0.17253734,  -0.11075424,
       -0.21586166,  0.10308852,   0.00525577,   0.09085813,   0.19607207,   0.25001204,
        0.03426333,  0.22383861,   0.04800471,   0.04972978,   0.09568697,   0.2430264,
        0.28215173,  0.17925027,   0.36056247,   0.33889654,   0.1248968 ,   0.16565846,
        0.23354977,  0.37686616,   0.2766894 ,   0.4223085 ,   0.57731646,   0.32448497,
        0.3310848 ,  0.48928034,   0.5434665 ,   0.5441702 ,   0.46737072,   0.60686255,
        0.5725527 ,  0.5917013 ,   0.46564984,   0.5967817 ,   0.68037146,   0.53769773,
        0.5776071 ,  0.6235754 ,   0.67287165,   0.82511437,   0.93292713,   0.84517425,
        0.762654  ,  0.8910726 ,   0.9227938 ,   0.7618559 ,   0.80128413]

ts1, ts2 = toNumpy(ts1, ts2)

# An example taken from P. Lepage's lsqfit tutorial:
ycov  = np.array(
       [[ 2.1537808808e-09,   8.8161794696e-10,   3.6237356558e-10,
          1.4921344875e-10,   6.1492842463e-11,   2.5353714617e-11,
          4.3137593878e-12,   7.3465498888e-13],
       [  8.8161794696e-10,   3.6193461816e-10,   1.4921610813e-10,
          6.1633547703e-11,   2.5481570082e-11,   1.0540958082e-11,
          1.8059692534e-12,   3.0985581496e-13],
       [  3.6237356558e-10,   1.4921610813e-10,   6.1710468826e-11,
          2.5572230776e-11,   1.0608148954e-11,   4.4036448945e-12,
          7.6008881270e-13,   1.3146405310e-13],
       [  1.4921344875e-10,   6.1633547703e-11,   2.5572230776e-11,
          1.0632830128e-11,   4.4264622187e-12,   1.8443245513e-12,
          3.2087725578e-13,   5.5986403288e-14],
       [  6.1492842463e-11,   2.5481570082e-11,   1.0608148954e-11,
          4.4264622187e-12,   1.8496194125e-12,   7.7369196122e-13,
          1.3576009069e-13,   2.3914810594e-14],
       [  2.5353714617e-11,   1.0540958082e-11,   4.4036448945e-12,
          1.8443245513e-12,   7.7369196122e-13,   3.2498644263e-13,
          5.7551104112e-14,   1.0244738582e-14],
       [  4.3137593878e-12,   1.8059692534e-12,   7.6008881270e-13,
          3.2087725578e-13,   1.3576009069e-13,   5.7551104112e-14,
          1.0403917951e-14,   1.8976295583e-15],
       [  7.3465498888e-13,   3.0985581496e-13,   1.3146405310e-13,
          5.5986403288e-14,   2.3914810594e-14,   1.0244738582e-14,
          1.8976295583e-15,   3.5672355835e-16]]
       )


mat = [[ 1,2,-1],
       [ 2,1, 2],
       [-1,2, 1]]
mat = np.array(mat)

psd =np.array(
[[ 1.        ,  0.64069129, -0.17902934],
 [ 0.64069129,  1.        ,  0.64069129],
 [-0.17902934,  0.64069129,  1.        ]])


def norm_cov(cov) -> np.ndarray:
    """ A more transparent covariance normalizer. """ 
    res = np.zeros((len(cov), len(cov[0])))
    for i in range(len(cov)):
        for j in range(len(cov[0])):
            res[i][j] = cov[i][j] / np.sqrt( cov[j][j] * cov[i][i] )
    return np.array(res)


def testStats():

    lpass = True

    A = []
    B = []
    for i in range(3):
        A.append(i)
        B.append(2-i)
    A = np.array(A)
    B = np.array(B)

    q12=gaudif(x1,e1,x2,e2)
    q13=gaudif(x1,e1,x3,e3)
    q14=gaudif(x1,e1,x4,e4)

    # Test gaussian difference
    if abs(q12-q12control) > eps:
        lpass=False
    if abs(q13-q13control) > eps:
        lpass=False
    if abs(q14-q14control) > eps:
        lpass=False

    s12=studif(x1,e1,ndat1,x2,e2,ndat2)

    # Test student difference
    if abs(s12-s12control) > eps:
        lpass=False

    # Student and Gaussian difference tests should agree for large sample sizes
    if studif(x1,e1,300,x2,e2,300)/gaudif(x1,e1,x2,e2) > 1.005:
        lpass=False

    # Explicit check that we have the correct interpretation of np.corrcoeff
    xbar1 = std_mean(ts1)
    xbar2 = std_mean(ts2)
    R = np.sum((ts1-xbar1)*(ts2-xbar2))/np.sqrt(np.sum((ts1-xbar1)**2)*np.sum((ts2-xbar2)**2))

    # Two simple tests of method normalizing a covariance matrix.
    lpass *= print_results(cov_to_cor(ycov),cov_to_cor(cov_to_cor(ycov)),text='cov_to_cor fixed point')
    lpass *= print_results(cov_to_cor(ycov),norm_cov(ycov),text='cov_to_cor')

    lpass *= print_results(pearson(ts1,ts2),R,text='pearson explicit')

    lpass *= print_results(psd,forcePositiveSemidefinite(mat),text='positive semidefinite',prec=1e-7)
    if not isPositiveSemidefinite(psd):
        logger.TBFail('psd not positive semidefinite')
        lpass = False

    # This will work whatever ddof is, since this factor cancels in the ratio.
    lpass *= print_results(cov_to_cor(np.cov(ts1,ts2))[0,1], pearson(ts1,ts2), text='pearson numpy')

    # Test of the confidence ellipse method.
    Ndraws = 1001
    rng = np.random.default_rng(seed=DEFAULTSEED)
    _,ax = plt.subplots()
    x = rng.normal(0,3,Ndraws)
    y = rng.normal(0,1,Ndraws)
    a,b,_ = confidence_ellipse(x,y,ax=ax,CI=0.5)
    Ninside=0
    for i in range(len(x)):
        if x[i]**2/a**2 + y[i]**2/b**2 < 1:
            Ninside+=1
    lpass *= print_results(Ninside/Ndraws,0.5004995004995005,text='ellipse')

    data1 = rng.normal(0,1,Ndraws)
    data2 = rng.normal(0,1,Ndraws)
    lpass *= print_results( KSTest_2side(data1,data2), 0.7359984219245113, text='KS 2 side')

    normalCDF = sp.stats.norm(loc=0,scale=1).cdf

    lpass *= print_results(KSTest_1side(data1,normalCDF),0.4937205464565033, text='KS 1 side 1')
    lpass *= print_results(KSTest_1side(data2,normalCDF),0.656583347109468 , text='KS 1 side 2')

    lpass *= print_results(covariance(ts1,ts1),std_var(ts1), text='diagonal cov')

    m, e = symmetrizeError(lo=0.1,hi=0.3,central=0.2)
    lpass *= print_results(m,0.2,e,0.3,text='conservative symmetric')

    m, e = symmetrizeError(lo=0.0067,hi=0.015,central=0.4529,method='FLAG')
    lpass *= print_results(m,0.454975,e, 0.012924999999999999,text='FLAG symmetric')

    concludeTest(lpass)


if __name__ == '__main__':
    testStats()