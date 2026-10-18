# 
# testSpline.py                                                               
# 
# D. Clarke
# 
# Test some of the convenience wrappers for spline methods.
# 

import numpy as np
import itertools
from math import comb
from latqcdtools.math.spline import _even_knots, _random_knots, getSpline, bootSpline, getSplineErr, modelSpread, \
    Interpolating, Regression, KnotAverage
from scipy.interpolate import CubicSpline, BSpline, splrep, splev
from latqcdtools.base.plotting import plt, plot_dots, plot_lines, set_params
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.statistics.statistics import countParams
import latqcdtools.base.logger as logger
from latqcdtools.base.initialize import DEFAULTSEED, TBRNG


SHOWPLOT = False


def testSpline():

    lpass = True

    x  = np.linspace(-1, 1, 101)
    y  = 10*x**2 + TBRNG(DEFAULTSEED).normal(size=len(x))
    ye = np.repeat(1.,len(y))

    knots = _even_knots(x, 3)

    lpass *= print_results(knots,[-0.5, 0.0, 0.5], text="even_knots")

    knots = _random_knots(x, 3, SEED=DEFAULTSEED)

    lpass *= print_results(knots,[-0.55, -0.10999999999999999, 0.10999999999999999], text="random_knots")

    # countParams must return the number of free spline parameters, not the padded
    # length of splrep's coefficient array (which overcounts by order+1) nor len(c)
    # for a scipy CubicSpline (which is just order+1).
    for nk in [2, 4, 6]:
        bspl = getSpline(x, y, ye, Regression(nknots=nk, order=3))
        # len(interior knots) + order + 1 free B-spline coefficients
        lpass *= print_results(countParams(bspl, ()), nk + 4, text=f"countParams B-spline nk={nk}")
    cspl = getSpline(x, y)
    lpass *= print_results(countParams(cspl, ()), len(x), text="countParams interpolating spline")

    # Natural least-squares spline: zero curvature at the endpoints, 2 fewer parameters.
    nspl = getSpline(x, y, ye, Regression(nknots=4, natural=True))
    lpass *= print_results(np.array([nspl(x[0],2),nspl(x[-1],2)]), np.zeros(2), text="natural LSQ S''=0 at ends", abs_prec=1e-10)
    lpass *= print_results(countParams(nspl, ()), 4 + 2, text="countParams natural LSQ spline")

    # With a knot at every interior datum, it must reproduce the natural interpolating spline.
    xs  = np.linspace(0, 1, 12)
    ys  = np.sin(3*xs)
    ispl = getSpline(xs, ys, np.ones(len(xs)), Regression(knots=xs[1:-1], natural=True))
    xf  = np.linspace(0, 1, 50)
    lpass *= print_results(ispl(xf), CubicSpline(xs, ys, bc_type='natural')(xf), text="natural LSQ -> interpolating limit")

    aic_arr = []
    for knots in [10,30,60]:

        spline = getSpline(x,y,ye,Regression(nknots=knots))
        aic    = spline.get_AIC()

        aic_arr.append(aic)

        if SHOWPLOT:
            set_params(xlabel='x',ylabel='y')
            plot_dots(x, y, ye)
            plot_lines(x,spline(x),marker=None)
            plt.show()

        logger.info('Try spline with knots =',knots,'AIC =',aic)

    if aic_arr[0]>aic_arr[-1]:
        logger.TBFail("AIC should penalize overfitting.")
        lpass=False
    else:
        logger.TBPass("AIC penalizes overfitting.")

    # bootSpline: the rebuilt bootstrap splines must equal refits to the same samples. For fixed knots the
    # fit is linear in y, so the band can be compared with exact error propagation.
    xb  = np.linspace(0, 3, 40)
    eb  = np.repeat(0.05,len(xb))
    yb  = np.sin(xb) + TBRNG(DEFAULTSEED).normal(0,0.05,len(xb))
    kind = Regression(nknots=5, natural=True)
    res = bootSpline(xb, yb, eb, kind, numb_samples=1000)
    sBS = res['splineBS'][0]
    lpass *= print_results(sBS(res['xspl']), getSpline(xb, sBS.get_y(), eb, kind)(res['xspl']),
                           text="bootSpline rebuilt spline = refit")
    lpass *= print_results(res['AIC'], getSpline(xb, yb, eb, kind).get_AIC(),
                           text="bootSpline AIC of original data")
    M = []
    for unit in np.eye(len(xb)):
        M.append(getSpline(xb, unit, eb, kind)(res['xspl']))
    M     = np.array(M).T
    exact = np.sqrt(np.diag(M @ np.diag(eb**2) @ M.T))
    lpass *= print_results(res['ysple'], exact, text="bootSpline error band", prec=0.15)
    lpass *= print_results(res['ymaxe'], res['ymaxe_stat'], text="fixed knots: no systematic error", prec=1e-12)
    def ymax(spl):
        return np.max(spl(res['xspl']))
    lpass *= print_results(modelSpread(res, ymax), 0., text="fixed knots: modelSpread is 0", abs_prec=1e-12)
    ymaxBS = []
    for spl in res['splineBS']:
        ymaxBS.append(ymax(spl))
    lpass *= print_results(res['ymax'], 2*ymax(res['splineMean']) - np.median(ymaxBS), text="ymax bias corrected", prec=1e-8)

    # Same with correlated data, which are now also fit with the covariance matrix: the band must follow 
    # sqrt(diag(M cov M^T)) with M the map from the data to the correlated fit.
    dist  = np.abs(np.subtract.outer(xb,xb))
    cov   = np.outer(eb,eb)*np.exp(-dist/0.5)
    res   = bootSpline(xb, yb, cov, kind, numb_samples=1000)
    M = []
    for unit in np.eye(len(xb)):
        M.append(getSpline(xb, unit, cov, kind)(res['xspl']))
    M     = np.array(M).T
    exact = np.sqrt(np.diag(M @ cov @ M.T))
    lpass *= print_results(res['ysple'], exact, text="bootSpline correlated error band", prec=0.15)

    # Regression: must agree with splrep (task=-1), with and without weights.
    xf = np.linspace(xb[0], xb[-1], 200)
    for nk in [1, 3, 7]:
        knots = _even_knots(xb, nk)
        lpass *= print_results(getSpline(xb, yb, eb, Regression(nknots=nk))(xf), splev(xf, splrep(xb, yb, w=1/eb, t=knots, k=3, task=-1)),
                               text=f"Regression = splrep nknots={nk}", prec=1e-8)
    lpass *= print_results(getSpline(xb, yb, kind=Regression(nknots=3, order=2))(xf), splev(xf, splrep(xb, yb, t=_even_knots(xb, 3), k=2, task=-1)),
                           text="unweighted quadratic Regression = splrep", prec=1e-8)

    # Correlated regression: compare with the normal equations c = (B^T W B)^-1 B^T W y.
    t     = np.r_[[xb[0]]*4, _even_knots(xb, 5), [xb[-1]]*4]
    Bmat  = BSpline.design_matrix(xb, t, 3).toarray()
    Winv  = np.linalg.inv(cov)
    c     = np.linalg.solve(Bmat.T @ Winv @ Bmat, Bmat.T @ Winv @ yb)
    lpass *= print_results(getSpline(xb, yb, cov, Regression(nknots=5))(xf), BSpline(t, c, 3)(xf), text="correlated Regression", prec=1e-8)

    # Interpolating: the natural cubic interpolant. It is the default kind, with and without edata.
    lpass *= print_results(getSpline(xb, yb, kind=Interpolating())(xf), CubicSpline(xb, yb, bc_type='natural')(xf), text="Interpolating = CubicSpline")
    lpass *= print_results(getSpline(xb, yb)(xf), CubicSpline(xb, yb, bc_type='natural')(xf), text="default kind without edata")
    lpass *= print_results(getSpline(xb, yb, eb)(xf), CubicSpline(xb, yb, bc_type='natural')(xf), text="default kind with edata")
    lpass *= print_results(bootSpline(xb, yb, eb, numb_samples=10)['splineMean'](xf), CubicSpline(xb, yb, bc_type='natural')(xf), text="bootSpline default kind")

    # getSplineErr with random knots: both fits must use the same knots, so the center is the fit to y.
    kind   = Regression(nknots=4, rand=True)
    center = getSplineErr(xb, xf, yb, eb, kind=kind)[0]
    lpass *= print_results(center, getSpline(xb, yb, kind=kind)(xf), text="getSplineErr same random knots", prec=1e-8)

    # KnotAverage: on a tiny problem, enumerate every configuration on the grid and compare the exact AIC-weighted
    # P(k) and average with the chain. This checks the birth/death proposal ratio.
    xt = np.linspace(0, 1, 9)
    et = np.full(9, 0.05)
    yt = np.sin(4*xt) + et*TBRNG(3).normal(size=9)
    xf = np.linspace(0, 1, 41)
    for natural in [False, True]:
        f  = KnotAverage(ngrid=8, natural=natural, nupdates=30000, nburn=1000, every=5).fitter(xt, et)
        My = f.M @ yt
        pk, curve, Z = np.zeros(f.kmax+1), np.zeros(len(xf)), 0
        for k in range(f.kmax+1):
            for idx in itertools.combinations(range(f.ngrid), k):
                fit = f._system(np.array(idx, dtype=int), My)
                if fit is None:
                    continue
                with np.errstate(under='ignore'):
                    w = np.exp(-0.5*fit['AIC'] + 20) / comb(f.ngrid, k)
                pk[k] += w
                curve += w*(BSpline.design_matrix(xf, fit['t'], 3).toarray() @ fit['c'])
                Z += w
        det = f.details(yt, xf)
        lpass *= print_results(det['kHistogram']/np.sum(det['kHistogram']), pk/Z, text=f"KnotAverage P(k) = exact, natural={natural}", abs_prec=0.03)
        lpass *= print_results(f.fit(yt)(xf), curve/Z, text=f"KnotAverage curve = exact, natural={natural}", abs_prec=0.005)

    # The averaged spline, expressed with a knot at every grid point, must equal the count-weighted average of
    # Regression fits with the saved configurations.
    kind = KnotAverage(nupdates=2000, nburn=200)
    spl  = getSpline(xb, yb, eb, kind)
    det  = kind.fitter(xb, eb).details(yb, xf)
    ref  = np.zeros(len(xf))
    for i in range(len(det['counts'])):
        ref += det['counts'][i]*getSpline(xb, yb, eb, Regression(knots=det['configs'][i]))(xf)
    lpass *= print_results(spl(xf), ref/np.sum(det['counts']), text="KnotAverage = average of Regression fits", prec=1e-8)

    # bootSpline reruns the chain per sample, seeded by the sample, so it does not depend on nproc.
    kind = KnotAverage(nupdates=500, nburn=100)
    res1 = bootSpline(xb, yb, eb, kind, numb_samples=8)
    res2 = bootSpline(xb, yb, eb, kind, numb_samples=8, nproc=2)
    lpass *= print_results(res1['ysple'], res2['ysple'], text="KnotAverage bootSpline nproc")
    det  = kind.fitter(xb, eb).details(yb, res1['xspl'])
    w    = det['counts']/np.sum(det['counts'])
    ycfg = []
    for cfg in det['configs']:
        ycfg.append(getSpline(xb, yb, eb, Regression(knots=cfg))(res1['xspl']))
    ycfg = np.array(ycfg)
    lpass *= print_results(res1['ysple_sys'], np.sqrt(w @ (ycfg - w @ ycfg)**2), text="KnotAverage ysple_sys = spread of Regression fits", prec=1e-8)
    lpass *= print_results(res1['ysple']**2, res1['ysple_stat']**2 + res1['ysple_sys']**2, text="KnotAverage ysple = stat + sys", prec=1e-8)
    def curve(spl):
        return spl(res1['xspl'])
    def ymax(spl):
        return np.max(spl(res1['xspl']))
    lpass *= print_results(modelSpread(res1, curve), res1['ysple_sys'], text="modelSpread reproduces ysple_sys", prec=1e-8)
    lpass *= print_results(modelSpread(res1, ymax), res1['ymaxe_sys'], text="modelSpread reproduces ymaxe_sys", prec=1e-8)
    lpass *= print_results(res1['splineMean'](xf), getSpline(xb, yb, eb, KnotAverage(nupdates=500, nburn=100))(xf),
                           text="KnotAverage bootSpline splineMean = getSpline")

    # Old positional call getSpline(x, y, nknots) must fail with a hint.
    try:
        getSpline(xb, yb, 3)
        logger.TBFail("getSpline with scalar edata should fail.")
        lpass = False
    except Exception:
        logger.TBPass("getSpline with scalar edata fails.")

    # bootSpline samples from edata, so edata=None must fail.
    try:
        bootSpline(xb, yb, None)
        logger.TBFail("bootSpline without edata should fail.")
        lpass = False
    except Exception:
        logger.TBPass("bootSpline without edata fails.")

    concludeTest(lpass)


if __name__ == '__main__':
    testSpline()
