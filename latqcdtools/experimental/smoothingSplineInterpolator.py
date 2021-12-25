import numpy as np
import sys
from tqdm import tqdm
from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline, splev, splrep
from scipy.optimize import minimize
import latqcdtools.logger as logger


def guess(x, y, k, s, w=None):
    """Do an ordinary spline fit to provide knots"""
    return splrep(x, y, w, k=k, s=s)

def err(c, x, y, t, k, w=None):
    """The error function to minimize"""
    diff = y - splev(x, (t, c, k))
    if w is None:
        diff = np.einsum('...i,...i', diff, diff)
    else:
        diff = np.dot(diff*diff, w**2)
    return np.abs(diff)

def constraint_spline(x, y, x0, k=3, s=None, w=None):
    if x0 is None:
        print("Error! No constraint provided!")
        exit()
    t, c0, k = guess(x, y, k, s, w=w)
    con = {'type': 'eq',
           'fun': lambda c: splev(x0, (t, c, k), der=1),
           }
    opt = minimize(err, c0, (x, y, t, k, w), constraints=con)
    copt = opt.x
    return UnivariateSpline._from_tck((t, copt, k))


class SmoothingSplineInterpolator:
    def __init__(self, dataContainer, method, order, use_sigma=True, 
            xmin=None, xmax=None, fromSample=False, smoothingFactor=None, rawSmoothingFactor=None, 
            constrainPos=None, disable_progressBar=False):
        self.method = method
        self.order = order
        self.use_sigma = use_sigma
        self.xmin = xmin
        self.xmax = xmax
        self.fromSample = fromSample
        self.smoothingFactor = smoothingFactor
        self.rawSmoothingFactor = rawSmoothingFactor
        self.x0 = constrainPos
        self.disable_progressBar = disable_progressBar

        means_x = dataContainer.means_x
        stds_x = dataContainer.stds_x

        if self.xmin is not None:
            self.xmin = (self.xmin - means_x)/stds_x
        if self.xmax is not None:
            self.xmax = (self.xmax - means_x)/stds_x

        
        self.dataContainer = dataContainer
        self.train_x = dataContainer.train_x
        self.train_y = dataContainer.train_y
        self.sigma_y = dataContainer.sigma_y
    

    def fit(self, xdata_interpol):
        

        all_fits = []
        for order in self.order:
            print("Fitting order", order)
            
            if self.fromSample:
                x = self.train_x[0]
            else:
                x = self.train_x


            if self.fromSample:
                max_sample = self.train_x.shape[0]
                progressBar = tqdm(range(max_sample), 
                                    disable=self.disable_progressBar,
                                            file=sys.stdout)
                samples = []
                for i in progressBar:
                    X = self.train_x[i]
                    Y = self.train_y[i]
                    if self.use_sigma:
                        sigma = self.sigma_y[i]
                    
                    if self.xmin is not None:
                        xmin = self.xmin[i]
                    else:
                        xmin = X[0]
                    if self.xmax is not None:
                        xmax = self.xmax[i]
                    else:
                        xmax = X[-1]

                    ind = (X <= xmax) & (X >= xmin)
                    
                    if len(X[ind]) <= order:
                        logger.TBError("len(x) < order with xmin={xmin[0]} and xmax={xmax[0]}")

                    if self.rawSmoothingFactor is None:
                        m = len(X[ind])
                        s_max = m + (2*m)**0.5
                        s_min = m - (2*m)**0.5
                        if self.smoothingFactor is not None:
                            smoothing = s_min + self.smoothingFactor * (s_max - s_min)
                        else:
                            smoothing = s_min
                    else:
                        smoothing = self.rawSmoothingFactor

                    if self.method == "smoothingSpline":
                        if self.use_sigma:
                            spl = UnivariateSpline(X[ind], Y[ind], w=1/sigma[ind], k=order, 
                                s=smoothing, check_finite=True)
                        else:
                            spl = UnivariateSpline(X[ind], Y[ind], k=order, 
                                s=smoothing, check_finite=True)
                    elif self.method == "constraintSmoothingSpline":
                        if self.use_sigma:
                            spl = constraint_spline(X[ind], Y[ind], self.x0, w=1/sigma[ind], k=order,
                                    s=smoothing)
                        else:
                            spl = constraint_spline(X[ind], Y[ind], self.x0, k=order,
                                    s=smoothing)
                    else:
                        print("Error! Unknown method!")
                        exit()

                    y = spl(xdata_interpol)
                    y_err_fake = np.zeros(len(y))
                    data = list([xdata_interpol,y,y_err_fake])
                    samples.append(data)
                    
                samples = np.asarray(samples)
                all_fits.append(samples)
            else:
                X = self.train_x
                Y = self.train_y
                if self.use_sigma:
                    sigma = self.sigma_y

                if self.xmin:
                    xmin=self.xmin[0]
                else:
                    xmin=X[0]
                if self.xmax:
                    xmax=self.xmax[0]
                else:
                    xmax=X[-1]

                ind = (X <= xmax) & (X >= xmin)
                
                if len(X[ind]) <= order:
                    logger.TBError("len(x) < order with xmin={xmin[0]} and xmax={xmax[0]}")
                
                if self.rawSmoothingFactor is None:
                    m = len(X[ind])
                    s_max = m + (2*m)**0.5
                    s_min = m - (2*m)**0.5
                    if self.smoothingFactor is not None:
                        smoothing = s_min + self.smoothingFactor * (s_max - s_min)
                    else:
                        smoothing = s_min
                else:
                    smoothing = self.rawSmoothingFactor

                if self.method == "smoothingSpline":
                    if self.use_sigma:
                        spl = UnivariateSpline(X[ind], Y[ind], w=sigma[ind], k=order, 
                            s=smoothing, check_finite=True)
                    else:
                        spl = UnivariateSpline(X[ind], Y[ind], k=order, 
                            s=smoothing, check_finite=True)
                elif self.method == "constraintSmoothingSpline":
                    if self.use_sigma:
                        spl = constraint_spline(X[ind], Y[ind], self.x0, w=sigma[ind], k=order,
                                    s=smoothing)
                    else:
                        spl = constraint_spline(X[ind], Y[ind], self.x0, k=order,
                                    s=smoothing)
                else:
                    print("Error! Unknown method!")
                    exit()

                y = spl(xdata_interpol)
                
                # Since UnivariateSpline does not return sigma...
                # This should be changed by computing the error manually
                y_err_fake = np.zeros(len(y)) 

                all_fits.append(np.stack((xdata_interpol,y,y_err_fake)))


        all_fits = np.asarray(all_fits)
        
        if self.fromSample:
            all_fits = np.transpose(all_fits,(0,2,1,3))
  
        for i in range(len(all_fits)):
            x = all_fits[i,0]
            y = all_fits[i,1]
            sigmaFit = all_fits[i,2]
            
            x,y,sigmaFit = self.dataContainer.postProcess(x,y, sigma=sigmaFit)
            all_fits[i,0] = x
            all_fits[i,1] = y
            all_fits[i,2] = sigmaFit
        
        if len(all_fits.shape) == 4:
            all_fits = np.transpose(all_fits, (0,2,1,3))
            s = all_fits.shape
            all_fits = np.reshape(all_fits,(s[0]*s[1],s[2],s[3]))
        x_sample = all_fits[:,0]
        y_sample = all_fits[:,1]
        sigma_sample = all_fits[:,2]
        return x_sample, y_sample, sigma_sample, smoothing