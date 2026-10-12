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
from latqcdtools.math.spline import _even_knots, _random_knots, getSpline, bootSpline, posteriorSpline, getSplineErr, \
    Interpolating, Regression, Smoothing, KnotAverage
from scipy.interpolate import CubicSpline, make_smoothing_spline, BSpline, splrep, splev
from scipy.linalg import null_space
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

    # Smoothing spline: must agree with scipy, with and without weights. Several data at one x act like one
    # datum with the summed weight and the weighted mean y, which scipy needs since it refuses repeated x.
    rng = TBRNG(DEFAULTSEED)
    xs  = np.sort(rng.uniform(0, 5, 30))
    es  = rng.uniform(0.05, 0.2, 30)
    ys  = np.sin(2*xs)/xs + rng.normal(0, es)
    xf  = np.linspace(xs[0], xs[-1], 200)
    sspl = getSpline(xs, ys, es, Smoothing(lam=0.1))
    lpass *= print_results(sspl(xf), make_smoothing_spline(xs, ys, w=1/es**2, lam=0.1)(xf), text="smoothing spline = scipy", prec=1e-8)
    lpass *= print_results(getSpline(xs, ys, kind=Smoothing(lam=0.1))(xf), make_smoothing_spline(xs, ys, lam=0.1)(xf), text="unweighted smoothing spline = scipy", prec=1e-8)
    lpass *= print_results(np.array([sspl(xs[0],2),sspl(xs[-1],2)]), np.zeros(2), text="smoothing spline S''=0 at ends", abs_prec=1e-10)
    xd   = np.r_[xs, xs[5]]
    yd   = np.r_[ys, ys[5]+0.1]
    ed   = np.r_[es, 0.1]
    order = np.argsort(xd, kind='stable')
    xd, yd, ed = xd[order], yd[order], ed[order]
    wm   = 1/es**2
    ym   = np.copy(ys)
    wm[5] = 1/es[5]**2 + 1/0.1**2
    ym[5] = (ys[5]/es[5]**2 + (ys[5]+0.1)/0.1**2)/wm[5]
    lpass *= print_results(getSpline(xd, yd, ed, Smoothing(lam=0.1))(xf), make_smoothing_spline(xs, ym, w=wm, lam=0.1)(xf), 
                           text="smoothing spline with repeated x", prec=1e-8)
    # Effective number of parameters runs from ndata (interpolation) to 2 (weighted straight line).
    lpass *= print_results(getSpline(xs, ys, es, Smoothing(lam=1e-12)).get_nparams(), len(xs), text="peff lam->0", prec=1e-6)
    lspl = getSpline(xs, ys, es, Smoothing(lam=1e12))
    lpass *= print_results(lspl.get_nparams(), 2, text="peff lam->inf", prec=1e-6)
    lpass *= print_results(lspl(xf), np.polyval(np.polyfit(xs, ys, 1, w=1/es), xf), text="lam->inf straight line", prec=1e-6)

    # Correlated smoothing spline. Compare with the normal equations c = (B^T W B + lam Omega)^-1 B^T W y, 
    # computing Omega independently with Simpson's rule on each interval, which is exact for B_i'' B_j''.
    covs = np.outer(es,es)*np.exp(-np.abs(np.subtract.outer(xs,xs))/0.7)
    lpass *= print_results(getSpline(xs, ys, np.diag(es**2), Smoothing(lam=0.1))(xf), sspl(xf), text="smoothing spline diagonal covariance", prec=1e-8)
    t  = np.r_[[xs[0]]*4, xs[1:-1], [xs[-1]]*4]
    nB = len(t) - 4
    xq, wq = [], []
    for i in range(len(xs)-1):
        a, b = xs[i], xs[i+1]
        xq += [a, (a+b)/2, b]
        wq += [(b-a)/6, 4*(b-a)/6, (b-a)/6]
    D2    = BSpline(t, np.eye(nB), 3).derivative(2)(np.array(xq))
    Omega = D2.T @ (np.array(wq)[:,None]*D2)
    Bmat  = BSpline.design_matrix(xs, t, 3).toarray()
    Winv  = np.linalg.inv(covs)
    c     = np.linalg.solve(Bmat.T @ Winv @ Bmat + 0.1*Omega, Bmat.T @ Winv @ ys)
    cspl  = getSpline(xs, ys, covs, Smoothing(lam=0.1))
    lpass *= print_results(cspl(xf), BSpline(t, c, 3)(xf), text="correlated smoothing spline", prec=1e-8)
    H = []
    for unit in np.eye(len(xs)):
        H.append(getSpline(xs, unit, covs, Smoothing(lam=0.1))(xs))
    lpass *= print_results(cspl.get_nparams(), np.trace(np.array(H).T), text="peff = tr(H)", prec=1e-8)

    # bootSpline with a correlated smoothing spline: the band must follow exact error propagation.
    res = bootSpline(xs, ys, covs, Smoothing(lam=0.1), numb_samples=1000)
    M   = []
    for unit in np.eye(len(xs)):
        M.append(getSpline(xs, unit, covs, Smoothing(lam=0.1))(res['xspl']))
    M   = np.array(M).T
    lpass *= print_results(res['ysple'], np.sqrt(np.diag(M @ covs @ M.T)), text="bootSpline correlated smoothing spline band", prec=0.15)
    lpass *= print_results(res['splineBS'][0].get_nparams(), cspl.get_nparams(), text="bootSpline smoothing spline peff", prec=1e-8)

    # Log evidence, compared up to a constant with the textbook REML formula: c = Z a + (prior covariance 
    # pinv(Omega)/lam part), flat prior on a, with Z spanning the straight lines.
    Z    = null_space(Omega, rcond=1e-10)
    Opi  = np.linalg.pinv(Omega, rcond=1e-10, hermitian=True)
    X0   = Bmat @ Z
    ours, ref = [], []
    for lam in [1e-3, 0.1, 10]:
        K   = covs + Bmat @ Opi @ Bmat.T / lam
        Ki  = np.linalg.inv(K)
        A   = X0.T @ Ki @ X0
        P   = Ki - Ki @ X0 @ np.linalg.solve(A, X0.T @ Ki)
        with np.errstate(under='ignore'):
            ref.append(-0.5*ys @ P @ ys - 0.5*np.linalg.slogdet(K)[1] - 0.5*np.linalg.slogdet(A)[1])
        ours.append(getSpline(xs, ys, covs, Smoothing(lam=lam)).logEvidence)
    ours, ref = np.array(ours), np.array(ref)
    lpass *= print_results(ours - ours[0], ref - ref[0], text="log evidence = REML", abs_prec=1e-8)

    # Averaging over lam: rebuild the average for the original data and for one bootstrap sample from
    # getSpline, with weights proportional to the evidence.
    res = bootSpline(xs, ys, es, Smoothing(nlam=40), numb_samples=100)
    lpass *= print_results(np.sum(res['weights']), 1., text="lam average weights normalized")
    for spl in [res['splineMean'], res['splineBS'][0]]:
        logEv, curves = [], []
        for lam in res['lams']:
            lspl = getSpline(xs, spl.get_y(), es, Smoothing(lam=lam))
            logEv.append(lspl.logEvidence)
            curves.append(lspl(res['xspl']))
        with np.errstate(under='ignore'):
            w = np.exp(np.array(logEv) - np.max(logEv))
        w[w < np.exp(-70)] = 0
        w /= np.sum(w)
        lpass *= print_results(spl(res['xspl']), w @ np.array(curves), text="lam average = evidence-weighted splines", prec=1e-8)
    k     = np.argmax(res['weights'])
    kspl  = getSpline(xs, ys, es, Smoothing(lam=res['lams'][k]))
    lpass *= print_results(res['splinesLam'][k](res['xspl']), kspl(res['xspl']), text="splinesLam = getSpline", prec=1e-8)
    lpass *= print_results(res['splinesLam'][k].logEvidence, kspl.logEvidence, text="splinesLam logEvidence", prec=1e-8)
    lpass *= print_results(res['splinesLam'][k].get_nparams(), kspl.get_nparams(), text="splinesLam peff", prec=1e-8)
    ylam  = []
    for spl in res['splinesLam']:
        ylam.append(spl(res['xspl']))
    ylam  = np.array(ylam)
    lpass *= print_results(np.sqrt(res['weights'] @ (ylam - res['weights'] @ ylam)**2), res['ysple_lam'], text="splinesLam reproduce ysple_lam", prec=1e-8)

    # posteriorSpline: the posterior covariance at fixed lam is (B^T W B + lam Omega)^-1 (Omega from 
    # Simpson above), the mean is the evidence-weighted average, and the samples reproduce the exact 
    # posterior variance (law of total variance) within Monte Carlo accuracy.
    res  = posteriorSpline(xs, ys, es, Smoothing(nlam=40), numb_samples=4000)
    k    = np.argmax(res['weights'])
    lpass *= print_results(res['splinesLam'][k].covPost, np.linalg.inv(Bmat.T @ np.diag(1/es**2) @ Bmat + res['lams'][k]*Omega), 
                           text="posterior covariance", prec=1e-6)
    lpass *= print_results(res['yspl'], res['splineMean'](res['xspl']), text="posterior mean = splineMean", prec=1e-8)
    ypost = []
    for spl in res['splineSamples']:
        ypost.append(spl(res['xspl']))
    ypost = np.array(ypost)
    lpass *= print_results(np.std(ypost, axis=0), res['ysple'], text="posterior samples reproduce ysple", prec=0.08)
    lpass *= print_results(res['ysple']**2, res['ysple_stat']**2 + res['ysple_sys']**2, text="posterior total variance")
    lpass *= print_results(res['yspl'], bootSpline(xs, ys, es, Smoothing(nlam=40), numb_samples=10)['splineMean'](res['xspl']),
                           text="posterior mean = bootSpline splineMean", prec=1e-8)

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

    # Interpolating: the natural cubic interpolant. Default kinds: Interpolating without edata, the average over
    # lam with edata.
    lpass *= print_results(getSpline(xb, yb, kind=Interpolating())(xf), CubicSpline(xb, yb, bc_type='natural')(xf), text="Interpolating = CubicSpline")
    lpass *= print_results(getSpline(xb, yb)(xf), CubicSpline(xb, yb, bc_type='natural')(xf), text="default kind without edata")
    lpass *= print_results(getSpline(xb, yb, eb)(xf), bootSpline(xb, yb, eb, numb_samples=10)['splineMean'](xf), text="default kind with edata", prec=1e-8)

    # getSplineErr with random knots: both fits must use the same knots, so the center is the fit to y.
    kind   = Regression(nknots=4, rand=True)
    center = getSplineErr(xb, xf, yb, eb, kind=kind)[0]
    lpass *= print_results(center, getSpline(xb, yb, kind=kind)(xf), text="getSplineErr same random knots", prec=1e-8)

    # posteriorSpline at fixed lam: the mean is the smoothing spline, and the band is exact.
    res   = posteriorSpline(xb, yb, eb, Smoothing(lam=1e-3), numb_samples=10)
    lpass *= print_results(res['yspl'], getSpline(xb, yb, eb, Smoothing(lam=1e-3))(res['xspl']), text="posterior fixed lam mean", prec=1e-8)
    Bspl  = BSpline.design_matrix(res['xspl'], res['splinesLam'][0].get_knots(), 3).toarray()[:,:len(xb)+2]
    lpass *= print_results(res['ysple'], np.sqrt(np.diag(Bspl @ res['splinesLam'][0].covPost @ Bspl.T)), text="posterior fixed lam band", prec=1e-8)

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
    lpass *= print_results(res1['splineMean'](xf), getSpline(xb, yb, eb, KnotAverage(nupdates=500, nburn=100))(xf),
                           text="KnotAverage bootSpline splineMean = getSpline")

    # Old positional call getSpline(x, y, nknots) must fail with a hint.
    try:
        getSpline(xb, yb, 3)
        logger.TBFail("getSpline with scalar edata should fail.")
        lpass = False
    except Exception:
        logger.TBPass("getSpline with scalar edata fails.")

    concludeTest(lpass)


if __name__ == '__main__':
    testSpline()
