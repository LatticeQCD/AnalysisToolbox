import re, os
import numpy as np
from matplotlib.lines import Line2D as L2D
import latqcdtools.readin as rd
from latqcdtools.fitting import IllegalArgumentError
from latqcdtools.statistics import mean_and_std_dev, std_mean, std_median, dev_by_dist
from latqcdtools.spline_interpolate import constraint_spline, multi_spline_fit
from latqcdtools.readin import read_in, read_in_pure
import latqcdtools.bootstr as bs
# Needed for plt
from latqcdtools.plotting import *


def linear(x, a, b):
    return a*x + b


class Extrapolator:

    _allowed_keys = ['xmin', 'xmax', 'plot_xmin', 'plot_xmax', 'no_tex', 'show_plot',
            'title', 'ylabel', 'xlabel', 'outname', 'plot_results', 'tol',
            'nsamples', 'nblocks', 'randomization_factor',
            'constraints', 'order', 'knots', 'nknots', 'base_point', 'extr_xdata', 'folder',
            'save_sample']

    def __init__(self, Nts, xdata, ydata, edata = None, data_input = None,
            method = 'gauss_btstr',
            **kwargs):

        self._init_args(kwargs)

        if data_input is None:
            if method == "btstr":
                data_input = 'raw'
            if method == "gauss_btstr":
                data_input = 'direct'
            if method == "from_sample":
                data_input = 'sample'
            if method == "direct":
                raise ValueError("Need data input (--data-input)")


        self._Nts = np.asarray(Nts)

        self._xdata = []

        # Averaged data
        self._av_ydata = []
        self._av_edata = []

        # Arrays of pure data on which bootstrap is run
        self._ydata = []
        self._edata = []

        self._method = method

        if data_input == 'direct' and (method == 'from_sample' or method == 'btstr'):
            if len(self._plot_results) == 0:
                raise ValueError("Direct data input does not work with from_sample or"
                        " btstr as method")

        if data_input == 'sample' and method == 'btstr':
            raise ValueError("Bootstrapping on a bootstrap sample is not a good idea")


        for i in range(len(xdata)):
            ind = (xdata[i] >= self._xmin) & (xdata[i] <= self._xmax)
            if any(ind):
                self._xdata.append(xdata[i][ind])
                self._ydata.append(ydata[i][ind])
                if edata is None:
                    self._edata.append(np.ones_like(self._ydata[i]))
                    self._no_edata = True
                else:
                    self._edata.append(edata[i][ind])
                    self._no_edata = False


        if self._nsamples is None:
            if data_input == "sample":
                    self._nsamples = len(self._ydata[i][0])
            else:
                self._nsamples = 1000


        if data_input == 'direct':
            self._av_ydata = self._ydata
            self._av_edata = self._edata

        elif data_input == 'sample':
            for yd in self._ydata:
                av_ydata, av_edata = mean_and_std_dev(yd, axis = 1)
                self._av_ydata.append(av_ydata)
                self._av_edata.append(av_edata)


        elif data_input == 'raw':
            for yd in self._ydata:
                av_ydata, av_edata = bs.bootstr(std_mean, yd, self._nsamples,
                        args = {'axis' : 1})

                self._av_ydata.append(av_ydata)
                self._av_edata.append(av_edata)


        self._constraints = np.asarray(kwargs.get('constraints', []))
        
        self._order = kwargs.get('order', 3)

        # Self defined knots
        self._knots = kwargs.get('knots', None)

        if self._knots is None:
            # Array of number of knots. For each number of knots a separate extrapolation is performed
            nknots = kwargs.get('nknots', None)
        else:
            nknots = [len(self._knots)]

        if nknots is None:
            if len(self._plot_results) > 1:
                self._nknots = []
            else:
                self._nknots = [4, 5, 6]
        else:
            self._nknots = nknots


        self._direct_res = []
        self._direct_res_err = []

        self._res_sample = []
        self._extr_sample = []

        self._direct_knots = []
        self._extr_data = None
        self._extr_data_err = None

        self._auto_plot_range()



    def _init_args(self, kwargs):
        diff = set(set(kwargs.keys()) - set(self._allowed_keys))
        if len(diff) != 0:
            raise IllegalArgumentError("Illegal argument(s) to fitter", *diff)


        self._xmin = kwargs.get('xmin', -np.inf)
        self._xmax = kwargs.get('xmax', np.inf)

        self._plot_xmin = kwargs.get('plot_xmin', None)
        self._plot_xmax = kwargs.get('plot_xmax', None)
        self._no_tex = kwargs.get('no_tex', False)

        self._show_plot = kwargs.get('show_plot', False)

        self._title = kwargs.get('title', None)
        self._ylabel = kwargs.get('ylabel', None)
        self._xlabel = kwargs.get('xlabel', None)
        self._folder = kwargs.get('folder', ".") + "/"
        os.makedirs(self._folder, exist_ok=True)
        self._outname = kwargs.get('outname', "extr")
        self._plot_results = kwargs.get('plot_results', [])

        self._base_point = kwargs.get('base_point', 0)

        self._tol = kwargs.get('tol', 1e-12)
        self._save_sample = kwargs.get('save_sample', False)

        self._nsamples = kwargs.get('nsamples', None)
        self._nblocks= kwargs.get('nblocks', 20)

        self._extr_xdata =  kwargs.get('extr_xdata', None)

        self._randomization_factor = kwargs.get('randomization_factor', 0.5)



    def _auto_plot_range(self):
        if self._plot_xmin is None:
            self._plot_xmin = np.inf
            for xd in self._xdata:
                if self._plot_xmin > np.min(xd):
                    self._plot_xmin = np.min(xd)

        if self._plot_xmax is None:
            self._plot_xmax = -np.inf
            for xd in self._xdata:
                if self._plot_xmax < np.max(xd):
                    self._plot_xmax = np.max(xd)

        if self._extr_xdata is None:
            self._extr_xdata = np.linspace(self._plot_xmin, self._plot_xmax, 1000)



    def perform_fits(self):
        self._res_sample = []
        self._direct_res = []
        self._direct_res_err = []


        for nknots in self._nknots:
            res, res_err, sample = self._compute_samples(nknots)

            self._res_sample.append(sample)
            self._direct_res.append(res)
            self._direct_res_err.append(res_err)


    def _compute_err_band(self, inv_Nt_squared):
        extr_sample = []
        for i in range(len(self._res_sample[0])):
            for j in range(len(self._res_sample)):
                res, knots = self._res_sample[j][i]
                coeffs = [linear(inv_Nt_squared, *j) for j in res]
                extr_data = constraint_spline(self._extr_xdata, knots, coeffs,
                        self._order, None, 0, self._constraints,
                        base_point = self._base_point)
                extr_sample.append(extr_data)


        extr_data = std_median(extr_sample, axis = 0)
        extr_data_err = dev_by_dist(extr_sample, axis = 0)
        return extr_sample, extr_data, extr_data_err

    def compute_err_band(self, inv_Nt_squared):
        extr_sample, extr_data, extr_data_err =  self._compute_err_band(inv_Nt_squared)
        return self._extr_xdata, extr_data, extr_data_err

    def compute_extrapolation(self, extr_xdata = None):

        if extr_xdata is not None:
            self._extr_xdata = extr_xdata

        self._extr_sample = []

        if not self._method == 'direct':
            self._extr_sample, self._extr_data, self._extr_data_err = self._compute_err_band(0)
        else:

            for i,res in enumerate(self._direct_res):
                coeffs = [linear(0, *j) for j in res]
                extr_data = constraint_spline(self._extr_xdata, self._direct_knots[i], coeffs,
                        self._order, None, 0, self._constraints, base_point = self._base_point)
                self._extr_sample.append(extr_data)

            self._extr_data = std_mean(self._extr_sample, axis = 0)

            # No error for direct method
            self._extr_data_err = np.full_like(self._extr_data, np.nan)

        return self._extr_xdata, self._extr_data, self._extr_data_err


    def get_extr_sample(self):
        return self._extr_sample




    def _fit_once(self, ydata, edata, start_params, nknots):
        knots, res, res_err, chi_dof = multi_spline_fit(linear, 2, 1/self._Nts**2,
                self._xdata, ydata, edata, order = self._order, tol = self._tol,
                always_return = True, base_point = self._base_point,
                constraints = self._constraints, start_params = start_params,
                nknots = nknots, algorithms = ["curve_fit"],
                randomization_factor = self._randomization_factor,
                knots = self._knots)


        return res, knots




    
    """
    To be called from bootstrap
    """

    def _fit_func(self, ydata, *args, **kwargs):
        current_ydata = []
        current_edata = []
        for i in range(len(ydata)):
            current_ydata.append(std_median(ydata[i], axis = 1))
            current_edata.append(dev_by_dist(ydata[i], axis = 1))
        return self._fit_once(current_ydata, current_edata, *args, **kwargs)



    def _compute_samples(self, nknots):

        # For start parameter estimation perform a direct fit
        knots, res, res_err, chi_dof = multi_spline_fit(linear, 2, 1/self._Nts**2, self._xdata,
                self._av_ydata, self._av_edata, randomization_factor = 0,
                order = self._order, tol = self._tol, always_return = True,
                constraints = self._constraints, base_point = self._base_point,
                nknots = nknots, knots = self._knots,
                algorithms = ["curve_fit"])

        logger.info("Reference knots of full fit = ", knots)
        self._direct_knots.append(knots)
        logger.info("Chi^2/d.o.f. of first full fit:", chi_dof)

        if self._method == 'gauss_btstr':
            sample, tmp_res, tmp_res_err = bs.bootstr_from_gauss(
                    self._fit_once, self._av_ydata, self._av_edata,
                    self._nsamples, err_by_dist = True, return_sample = True, 
                    args = {'edata' : self._av_edata, 'start_params' : res, 'nknots' : nknots},
                    nmax_exceptions = 0.25 * self._nsamples)

        elif self._method == 'from_sample':
            sample = []
            
            nsamples = len(self._ydata[0][0])
            nsteps = min(self._nsamples, nsamples)

            if nsteps < nsamples:
                logger.warn("Not using all samples that are available.")

            if nsteps < 100:
                step = 1
            else:
                step = nsteps / 100

            for i in range(nsteps):

                current_ydata = []
                current_edata = []
                for j in range(len(self._ydata)):
                    current_ydata.append(self._ydata[j][:,i])
                    if self._no_edata:
                        current_edata.append(self._av_edata[j])
                    else:
                        current_edata.append(self._edata[j][:,i])

                tmp_res = self._fit_once(current_ydata, current_edata,
                        res, nknots)
                sample.append(tmp_res)
                if i%step == 0:
                    logger.progress("%d%%" % ((i+1)/nsteps * 100))
        
        elif self._method == 'btstr':
                
            sample, tmp_res, tmp_res_err = bs.bootstr(
                    self._fit_func, self._ydata, self._nsamples,
                    conf_axis = 2, same_rand_for_obs = True,
                    err_by_dist = True, return_sample = True, 
                    args = {'start_params' : res, 'nknots' : nknots},
                    nmax_exceptions = 0.25 * self._nsamples)

        elif self._method == 'direct':
            sample = None

        else:
            raise ValueError("No such method: " + str(self._method))


        return res, res_err, sample



    def plot_extrapolation(self, outname = None, plot_size = (15, 10), kwargs_fill = {}, kwargs_dots = {}):

        if outname is not None:
            if self._no_tex:
                init_notex(fig_width = plot_size[0], fig_height = plot_size[1])
            else:
                latexify(fig_width = plot_size[0], fig_height = plot_size[1])

        kwargs_fill.setdefault('alpha', 0.5)

        line_styles = ['-', ':', '--', '-.']

        for i in range(len(self._xdata)):
            plot_dots(self._xdata[i], self._av_ydata[i], self._av_edata[i],
                    label = "$N_\\tau$ = %d" %(self._Nts[i],),
                    color = colors[i], xmin = self._plot_xmin, xmax = self._plot_xmax, **kwargs_dots)


        for i, Nt in enumerate(self._Nts):
            try:
                sample, ydata, edata = self._compute_err_band(1/Nt**2)
                plot_fill(self._extr_xdata, ydata, edata,
                        xmin = self._plot_xmin, xmax = self._plot_xmax,
                        color = colors[i], **kwargs_fill)
                

            except (IndexError, TypeError):
                # If no sample is available
                for j, res in enumerate(self._direct_res):
                    coeffs = [linear(1/self._Nts[i]**2, *j) for j in res]

                    plot_func(constraint_spline, xmin = self._plot_xmin, xmax = self._plot_xmax,
                            args = (self._direct_knots[j], coeffs, self._order, None, 0,
                                self._constraints, self._base_point),
                            color = colors[i], alpha = 0.0,
                            linestyle = line_styles[j%len(line_styles)])



        if self._extr_data is not None:
            plot_fill(self._extr_xdata, self._extr_data, self._extr_data_err,
                    xmin = self._plot_xmin, xmax = self._plot_xmax,
                    color = 'black', label = "cont.", **kwargs_fill)
        else:
            for i in range(len(self._direct_res)):
                handles, labels = get_legend_handles()
                handles.append(L2D([1, 2], [1, 2], linestyle = line_styles[i % len(line_styles)],
                    linewidth = 0.5, color = 'gray'))
                labels.append("nknots = " + str(self._nknots[i]))
            for i, res in enumerate(self._direct_res):
                coeffs = [linear(0, *i) for i in res]
                plot_func(constraint_spline, xmin = self._plot_xmin, xmax = self._plot_xmax,
                        args = (self._direct_knots[i], coeffs, self._order,
                            None, 0, self._constraints, self._base_point),
                        color = "black", alpha = 0.0,
                        linestyle = line_styles[i%len(line_styles)])



        
        set_params(xlabel = self._xlabel, ylabel = self._ylabel, title = self._title,
                legendpos = "best")

        if outname is not None:
            plt.savefig(self._folder + outname + ".pdf")

        if self._show_plot:
            plt.show()




    def save_extrapolation(self, outname = None):
        if outname is not None:
            self._outname = outname


        with open(self._folder + self._outname + "_parameters.txt", "w") as fout:
            print("#method:", self._method, file = fout)
            print("#base_point:", self._base_point, file = fout)
            print("#constraints:", file = fout)
            for cons in self._constraints:
                print(*cons, file = fout)
            for i, nknots in enumerate(self._nknots):
                print("#order = %d, nknots = %d, Tmin = %f, Tmax = %f"
                        % (self._order, nknots, self._xmin, self._xmax), file = fout)
                print("#parameters:", file = fout)
                for j in self._direct_res[i]:
                    print(j[0], j[1], file = fout)

                print("#knots:", file = fout)
                for k in self._direct_knots[i]:
                    print(k, file = fout)


        with open(self._folder + self._outname + "_cont.txt", "w") as fout:

            print("#xdata extr_ydata extr_ydata_err", file = fout)

            for j in range(len(self._extr_xdata)):
                print(self._extr_xdata[j], self._extr_data[j],
                        self._extr_data_err[j], file = fout)

        for i, Nt in enumerate(self._Nts):
            with open(self._folder + self._outname + "_%i.txt" % Nt, "w") as fout:
                sample, ydata, edata = self._compute_err_band(1/Nt**2)
                for j in range(len(self._extr_xdata)):
                    print(self._extr_xdata[j], ydata[j],
                            edata[j], file = fout)
                

        if self._method != "direct":
            with open(self._folder + self._outname + "_coeffs.txt", "w") as fout:
                print("#sample_index knots / coeffs", file = fout)
                for n in range(len(self._res_sample)):
                    print("#nknots = ", len(self._direct_knots[n]), file = fout)
                    for i in range(len(self._res_sample[n])):
                        print(i, *self._res_sample[n][i][0].flatten(), end = " / ", file = fout)
                        print(*self._res_sample[n][i][1], file = fout)

        if self._save_sample:
            with open(self._folder + self._outname + "_sample.txt", "w") as fout:
                print("#xdata extr_ydata", end = " ", file = fout)

                int_ydata = []
                if self._method != "direct":
                    for Nt in self._Nts:
                        print("#int_ydata_Nt" + str(Nt), end = " ", file = fout)
                        int_ydata.append(self._compute_err_band(1/Nt**2)[0])

                    fout.write("\n")

                for j in range(len(self._extr_sample)):
                    for i in range(len(self._extr_xdata)):
                        print(self._extr_xdata[i], self._extr_sample[j][i], end = " ",
                                file = fout)

                        if self._method != "direct":
                            for k in range(len(self._Nts)):
                                print(int_ydata[k][j][i], end = " ", file = fout)
                            fout.write("\n")



    def read_extrapolation(self, extr_xdata = None):

        append_to_params = False
        append_to_knots = False
        append_to_constraints = False
        self._constraints = []

        with open(self._plot_results[0]) as fin:
            for line in fin:
                if line.startswith("#parameters:"):
                    self._direct_res.append([])
                    append_to_params = True
                    append_to_knots = False
                    append_to_constraints = False
                    continue
                if line.startswith("#constraints:"):
                    append_to_params = False
                    append_to_constraints = True
                    append_to_knots = False
                    continue
                if line.startswith("#knots:"):
                    self._direct_knots.append([])
                    append_to_params = False
                    append_to_constraints = False
                    append_to_knots = True
                    continue
                if line.startswith('#order'):
                    self._order = int(line.split()[2].split(",")[0])
                    continue
                if line.startswith('#base_point:'):
                    self._base_point = float(line.split()[1])
                    continue

                if append_to_params:
                    self._direct_res[-1].append([float(i) for i in line.split()])
                if append_to_constraints:
                    self._constraints.append([float(i) for i in line.split()])
                if append_to_knots:
                    self._direct_knots[-1].append(float(line))

        self._nknots = [len(i) for i in self._direct_knots]

        if len(self._plot_results) > 1:
            self._extr_xdata, self._extr_data, self._extr_data_err = read_in(
                    self._plot_results[1], 1, 2, 3)

        if len(self._plot_results) > 2:
            self._res_sample = []
            with open(self._plot_results[2]) as fin:
                for line in fin:
                    if line.startswith("#nknots"):
                        self._res_sample.append([])

                    if not line.startswith('#'):
                        coeffs = [float(i) for i in line.split('/')[0].split()[1:]]
                        coeffs = np.array(coeffs).reshape(-1, 2)
                        knots = np.array([float(i) for i in line.split('/')[1].split()])
                        self._res_sample[-1].append((coeffs, knots))

            if extr_xdata is not None:
                self._extr_xdata = extr_xdata
            else:
                self._extr_xdata = np.linspace(self._plot_xmin, self._plot_xmax, 1000)
            self.compute_extrapolation()



def read_in_extr_files(files, xdata_col = 1, ydata_col = 2, edata_col = 3, Nts = None,
        read_method = 'raw'):
    if Nts is None:
        Nts = []
    new_Nts = []
    xdata = []
    ydata = []
    edata = []
    for filename in files:
        if read_method == 'direct':
            if edata_col is None:
                edata_col = 3

            tmp_x, tmp_y, tmp_e = rd.read_in(filename, xdata_col, ydata_col,
                    edata_col)
            if len(tmp_x) == 0:
                raise ValueError("File" + filename + "is empty")

            xdata.append(tmp_x)
            ydata.append(tmp_y)
            edata.append(tmp_e)

        elif read_method == 'sample' or read_method == 'raw':
            tmp_x, tmp_y, nconfs = read_in_pure(filename, xdata_col, ydata_col)
            xdata.append(tmp_x)
            ydata.append(tmp_y)
            if edata_col is None:
                edata.append(np.ones_like(tmp_y))
            else:
                tmp_x, tmp_e, nconfs = read_in_pure(filename, xdata_col, edata_col)
                edata.append(tmp_e)

        if len(Nts) == 0:
            fcheck=re.compile('(.*?_[Nn]t)(.*?)(\D+.*?)')
            x = fcheck.match(filename)
            if x is not None:
                Nt = x.group(2)
                new_Nts.append(int(Nt))
    
    if len(Nts) == 0:
        Nts = np.asarray(new_Nts)

    if len(Nts) == 0:
        raise ValueError("Could not get Nt from filename. Please pass with --Nts")

    if len(xdata) != len(Nts):
        raise ValueError("Number of filenames does not agree with number of Nts")
        

    return xdata, ydata, edata, Nts