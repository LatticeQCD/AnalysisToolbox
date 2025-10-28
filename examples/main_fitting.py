# 
# fitExample.py
# 
# D. Clarke, H. Sandmeyer
# 
# A basic fit example followed by a two-parameter fit example.
#


import numpy as np
from latqcdtools.base.plotting import plt, clearPlot, plot_dots, plot_lines, getColorGradient
from latqcdtools.statistics.fitting import Fitter, zipXYData, unzipXYData
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.initialize import initialize, finalize


initialize()


# Here we define our fit function. we pass it its independent variable followed by the fit parameters we are
# trying to determine.
def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*x**2 + b*x + c


xdata, ydata = readTable("../datasets/wurf.dat", usecols=(0,2)) 

# We initialize our Fitter object. The function has to look like
#            func(x, params, *args).
fitter = Fitter(fit_func, xdata, ydata)

# Here we try a fit, using the 'curve_fit' method, specifying the starting guesses for the fit parameters. Since
# detailedInfo = True, we will get back additional statistical information like the log of the Gaussian Bayes factor. 
res, res_err, chi_dof = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['curve_fit'], show_results=True)

# We can plot the fit and data using commands like this one. You can combine these commands with anything from
# plotting.py if you want to spice up your plot a bit.
fitter.plot_fit(domain=(np.min(xdata),np.max(xdata)))
fitter.plot_data()
plt.show()
clearPlot()



# OK, now let's try to fit a function of two variables.


# We start by making some fake points.
x  = np.linspace(-1,1,21)
y  = np.linspace(-1,1,5) 
Nx = len(x)
Ny = len(y)

# The zipXYData function conveniently prepares these x and y data
# for use inside an aribrary fit function of two variables.
xydata = zipXYData(x,y)

def fitFunc(data,c):
    # And here we can extract X and Y arrays. The zip/unzip lets us
    # pass data as one array, which is needed for the fit methods.
    X, Y = unzipXYData(data) 
    return c[0] + c[1]*X**2 + c[2]*Y

# We will have Z = 1 + 2x^2 + 3y. Add some fake noise.
rng = np.random.default_rng()
z = fitFunc(xydata,[1,2,3]) + rng.normal(0,0.1,len(xydata))

# Now we pass xydata as our array of independent data.
fitter = Fitter(fitFunc,xydata,z)
res, res_err, chi_dof = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['curve_fit'], show_results=True)

# Finally we plot the data. We decide to show cross-sections of f(x,y) for various fixed y.
Nplot = 101
colors = getColorGradient(Ny)
xplot  = np.linspace(-1,1,Nplot)
xyplot = zipXYData(xplot,y)
for i in range(Ny):
    # The data are ordered such that the first Nx data have the same y-value. Then the
    # next Nx data have a different y-value, and so on.
    zdata  = z     [i*Nx:(i+1)*Nx] 
    xdata  = xydata[i*Nx:(i+1)*Nx][:,0]
    ylabel = xydata[i*Nx:(i+1)*Nx][:,1][0]
    plot_dots(xdata,zdata,marker='o',color=colors[i],label=str(ylabel))
    xyfit  = xyplot[i*Nplot:(i+1)*Nplot]
    plot_lines(xplot,fitFunc(xyfit,res),marker=None,color=colors[i])

plt.show()
finalize()


