latqcdtools.math.spline
=============

`_even_knots(xdata, nknots)`
 
    Return a list of nknots evenly spaced knots. 
    
`_random_knots(xdata, nknots, randomization_factor=1, SEED=None)`
 
    Return a list of nknots randomly spaced knots. 
    
`getSpline(xdata, ydata, num_knots=None, edata=None, order=3, rand=False, fixedKnots=None, getAICc=False, natural=False, smooth=None)`
 
    This is a wrapper that calls SciPy spline fitting methods, depending on your needs. Calls LSQUnivariateSpline
    by default. If you need to ensure a well defined second derivative at the knots, we call instead UnivariateSpline,
    since LSQUnivariate spline seems to have no smoothing option. Sadly if you call UnivariateSpline, you can't specify
    the knots, so you can't both pick knots and smooth.

    Args:
        xdata (array-like)
        ydata (array-like)
        num_knots (int):
            The number of knots.
        edata (array-like, optional): 
            Error data. Defaults to None.
        order (int, optional):
            Order of the spline. Defaults to 3.
        rand (bool, optional): 
            Use randomly placed knots? Defaults to False.
        fixedKnots (array-like, optional):
            List of user-specified knots. Defaults to None.
        getAICc (bool, optional): 
            Return corrected Aikake information criterion? Defaults to False.
        natural (bool, optional): 
            Try a natural (no change in slope at the endpoints) cubic spline. Defaults to False. 

    Returns:
        callable spline object
        AICc (optionally)
    
`getSplineErr(xdata, xspline, ydata, ydatae, num_knots=None, order=3, rand=False, fixedKnots=None, natural=False)`
 
    Use getSpline to smooth mean and error bars. Create a spline-smooth band from that. 
    
`customSpline(xdata, ydata, edata=None, knots=None, order=3, smooth=None)`


