#!/usr/bin/env python3
import os, traceback, math
from latqcdtools.jackknife import jackknife
import latqcdtools.bootstr as bs
import latqcdtools.logger as logger
from latqcdtools.corr_fitting import CorrFitter
from latqcdtools.fitting import print_res, print_scl
from latqcdtools.corr_ratio_fitting import corr_ratio, corr_ratio_direct, CorrFitterRatio
from latqcdtools.readin import read_in, read_in_pure, read_in_fitmass
from latqcdtools.write import write_eff_mass, write_fit_mass, write_pcov, write_corr, write_sample
from latqcdtools.effmass import calc_eff_mass, calc_eff_mass_direct
from latqcdtools.remove_osc import remove_osc, remove_osc_av
try:
    from sklearn.covariance import MinCovDet as MCD
    sklearn_avail=True
except ImportError:
    sklearn_avail=False
import numpy as np
from latqcdtools.statistics import calc_cov, std_mean, std_dev, mean_and_err, mean_and_std_dev
from latqcdtools.tools import remove_nan


class ParameterError(Exception):
    pass


class PerformCorrFit:

    # Default start value for the lower value of the fit interval
    _nt_min_start = 1

    default_args = {
            # For detailed discription see bin/corrfitter.py --help
            'filename' : None, # Input file name

            'nstates' : 1, # Number non-oscillating states
            'nstates_osc' : 0, # Number oscillating states
            'sym' : True, # Symmetrize correlator directly after readin
            'change_sign' : False, # Multiply correlator with -(-1)**nt
            'auto_sign' : False, # Set above flag automatically
            'cut_eig' : False, # Cut eigenvalues of covariance matrix
            "cut_perc" : 30, # Percentage of eigenvalues to be cut
            'correlated' : False, # Perform a correlated fit
            'min_cov_det' : False, # Use minimal covariance determinat method
            # Support fraction of MCD method.
            # See http://scikit-learn.org/stable/modules/generated/sklearn.covariance.MinCovDet.html
            'mcd_supp_frac' : None,
            'numb_samples' : 100, # Number of samples for bootstrap for data operations
            'nfit_samples' : 100, # Number of samples for bootstrap for fits
            'numb_blocks' : 20, # Number of blocks for Jackknife

            'xmax' : None, # Upper Limit of fit interval
            'fit_interval' : None, # Range between nt_start is varied
            'priorval' : None, # Prior values for constraint fit
            'priorsigma' : None, # Prior sigma for constraint fit
            'seed' : None, # Seed for boostrap
            'try_all' : None, # Try all different start parameters. Expensive
            'start_params' : None, # Star parameters for fit
            'plot_start' : False, # Plot start parameters

            'folder' : "results", # Output folder
            'file_string' : "", # File string to be added to the output files
            'nt_col' : 1, # Column for nt in input file
            'data_col' : 2, # Column for data in input file
            'err_col' : 3, # For direct readin: Column for errors in input file

            'jack_data' : False, # Use Jackknife correlator data
            'std_err_data' : False, # Use standard error for correlator data
            'sample_data' : False, # Correlator data stem from bootstrap sample
            'direct_data' : False, # Read in already averaged data
            'btstr_data' : False, # Use bootstrapped correlator data
            'ratio_data' : False, # Use G(nt+1)/G(nt) as data

            'jack_fit' : False, # Estimate error using Jackknife
            'direct_fit' : False, # Estimate error directly from fit
            'btstr_fit' : False, # Estimate error using booststrap
            'ratio_fit' : False, # Fit the ration G(nt+1)/G(nt) instead of G(nt)
            'ng_btstr' : True, # Use quantiles for error distribution (non-gaussian -> ng)
            'scnd_btstr' : False, # Use second boostrap to get errors for fit in bootstrap fit

            'notex' : False, # Do not use Latex for plot rendering
            'xlabel' : "$n_{\\tau/\\sigma,\mathrm{min.}}$", # Default xlabel
            'ylabel' : "", # Default ylabel. If empty is set later on
            'title' : None, # Plot title
            'plot_no_ylog' : False, # Do not scale y-axis with log
            'plot_size' : (18,12), # Default plot size
            'font_size' : 11, # Default font size

            # Plot already performed fit results. This is the name of the input file
            'res_filename' : None,

            # Additionlly plot the correlator normalized with Gfree. This is m/T of Gfree
            'Nt': None # Number of correlator values.
            }


    def _check_args(self):
      # If we do not have a fit flag switch to default direct fit
        if not (self._jack_fit or self._btstr_fit or self._direct_fit or self._ratio_fit
                or self._plot_start or self._res_filename is not None):
            self._direct_fit = True


        # check if we did not get an argument for reading data
        if not (self._jack_data or self._btstr_data or self._direct_data or self._std_err_data or
                self._sample_data or self._ratio_data):

            # Generate read in flag from fit flag
            if self._jack_fit or self._plot_start:
                self._jack_data = True

            if self._btstr_fit:
                self._btstr_data = True

            if self._ratio_fit:
                self._ratio_data = True

            if self._direct_fit :
                # the default correlated fit shoud use standard error data
                if self._correlated:
                    self._std_err_data = True
                else:
                    self._jack_data = True

        # If there is still no read_in flag, set it direct
        if not (self._jack_data or self._btstr_data or self._direct_data or self._std_err_data or
                self._sample_data or self._ratio_data):
            self._std_err_data = True

        if self._ratio_fit and any([self._jack_data, self._btstr_data, self._direct_data,
            self._std_err_data, self._sample_data]):
            raise ParameterError("Ratio fit can only be combined with ratio_data")

        if self._ratio_fit and self._nstates_osc > 0:
            raise ParameterError("Ratio fit does only work for --nstates-osc 0 ")

        if np.sum([self._jack_data, self._btstr_data, self._direct_data, self._std_err_data,
            self._sample_data, self._ratio_data]) > 1:
            raise ParameterError("More than one read-in flag")

        if np.sum([self._jack_fit, self._btstr_fit, self._direct_fit, self._ratio_fit]) > 1:
            raise ParameterError("More than one fit flag")

        if self._direct_data and not (self._direct_fit or self._res_filename is not None):
            raise ParameterError("Direct data does only work with direct fit")

        if self._direct_data and self._correlated:
            raise ParameterError("Cannot compute covariance matrix if we do not know the whole"
                " data set: Do not use --direct-data")

        if self._sample_data and (self._btstr_fit or self._jack_fit):
            raise ParameterError("Data from bootstrap sample cannot be combined with jackknife"
            " or (non-gauss) bootstrap fit")


        if self._correlated and self._scnd_btstr:
            raise ParameterError("Correlated fits not possible with second bootstrap")

        if self._seed and self._scnd_btstr:
            raise ParameterError("Using a fixed seed and a second bootstrap is not supported")


        if self._min_cov_det and self._ratio_fit:
            raise ParameterError("Ratio fit and MinCovDet is not supported at the same time")


        if self._priorsigma is None and self._priorval is not None:
            raise ValueError("--prior-sigma has to be passed along with --prior-val")

        if self._priorsigma is not None and self._priorval is None:
            raise ValueError("--prior-val has to be passed along with --prior-sigma")


        if not self._ratio_fit:
            if self._ylabel == "":
                self._ylabel = "$G(n_{\\tau/\\sigma})$"
        else:
            self._ylabel = "$G(n_{\\tau/\\sigma})/G(n_{\\tau/\\sigma} + 1)$"

        if self._plot_no_ylog or self._ratio_fit:
            self._plot_ylog = False
        else:
            self._plot_ylog = True





    def _init(self):
        if self._seed is not None:
            self._out_name += "_seed" + str(self._seed)

        if self._correlated:
            self._out_name += "_cov"

        if self._min_cov_det:
            self._out_name += "_mcd"

        if not self._ng_btstr:
            self._out_name += "_no_median"

        if self._jack_data:
            self.read_in = self.read_in_jk
            self._out_name += "_jk-data"

        elif self._btstr_data:
            self.read_in = self.read_in_bs
            self._out_name += "_bs-data"

        elif self._std_err_data:
            self.read_in = self.read_in_std_err
            self._out_name += "_std-err-data"

        elif self._sample_data:
            self.read_in = self.read_in_sample
            self._out_name += "_fr-sample"

        elif self._direct_data:
            self.read_in = self.read_in_direct
            self._out_name += "_direct-data"

        elif self._ratio_data:
            if self._correlated:
                self.read_in = self.read_in_std_err_ratio
            else:
                self.read_in = self.read_in_jk_ratio
            self._out_name += "_ratio-data"

        # Read in data and compute errors
        self.read_in()

        # Compute the effective mass and also split correlators for nstates_osc > 0
        self._compute_meff()

        if self._change_sign:
            self._out_name += "_sc"

        if self._cut_eig:
            self._out_name += "_cut_" + str(self._cut_perc)


        self.init_fitter()


        if self._direct_fit:
            self.fit = self.direct_fit
            self._out_name += "_direct-fit"
            if self._try_all is None:
                self._try_all = True

        elif self._btstr_fit:
            self.fit = self.bs_fit

            if self._scnd_btstr:
                self._out_name += "_scnd-bs-fit"
            else:
                self._out_name += "_bs-fit"

            if self._try_all is None:
                self._try_all = False



        elif self._jack_fit:
            self.fit = self.jk_fit
            self._out_name += "_jk-fit"

            if self._try_all is None:
                self._try_all = False


        if self._ratio_fit:
            self.fit = self.direct_fit
            self._out_name += "_ratio-fit"

            if self._try_all is None:
                self._try_all = False





    def __init__(self, **kwargs):

        for key in kwargs:
            if key not in self.default_args:
                raise ParameterError("Key " + key + " not a valid parameter")
        for key, value in self.default_args.items():
            if key in kwargs:
                setattr(self, "_" + key, kwargs[key])
            else:
                setattr(self, "_" + key, value)
        self._check_args()

        if self._min_cov_det:
            if not sklearn_avail:
                raise ImportError("To use minimal covariance determinant method,"
                        "you need to install sklearn")


        # Initialize variables
        self._init_arrays()
        self._nparams = 2 * (self._nstates + self._nstates_osc)

        if self._ratio_fit:
            self._nparams -= 1

        self._corr_fitter = None
        self._osc_corr_fitter = None
        self._no_corr_fitter = None

        if self._nstates_osc > 0:
            self._osc = True
        else:
            self._osc = False

        if self._start_params is not None:
            if len(self._start_params) != self._nparams:
                raise ValueError("Number of start parameters does"
                        "not match to number of states")

        if self._sym:
            self._out_name = "sym" + self._file_string
        else:
            self._out_name = "asym" + self._file_string


        os.makedirs(self._folder, exist_ok=True)
        np.seterr(invalid = 'ignore')
        np.seterr(over = 'ignore')

        self._init()


    def _init_arrays(self):
        self._data = []
        self._rat_data = []

        self._xdata = []
        self._ydata = []
        self._edata = []

        self._xdata_osc = []
        self._ydata_osc = []
        self._edata_osc = []

        self._xdata_no = []
        self._ydata_no = []
        self._edata_no = []

        self._cov = []
        self._ratio_cov = []

        self._res = []
        self._res_err = []
        self._pcov = []
        self._chi_dof = []
        self._aicc = []
        self._ranges = []

        self._sep = []
        self._sep_err = []
        self._chi_dof_sep = []
        self._aicc_sep = []

        self._res_osc = []
        self._res_osc_err = []
        self._chi_dof_osc = []
        self._aicc_osc = []

        self._res_no = []
        self._res_no_err = []
        self._chi_dof_no = []
        self._aicc_no = []

        self._samples = []

        self._meff = []
        self._meff_err = []
        self._meff_no = []
        self._meff_no_err = []
        self._meff_osc = []
        self._meff_osc_err = []




    def _calc_fit_interval(self):

        #  stopnumb indicates how many points are needed for the fit.
        # In case of an asymetric fit we have 2*stopnumb data points
        if self._sym:
            stopnumb = 2 * (self._nstates + self._nstates_osc)
        else:
            stopnumb = (self._nstates + self._nstates_osc)

        if self._fit_interval is None:
            # This is a good estimate for the best interval.
            # Feel free to change
            self._fit_interval = [self._nt_min_start, min(self.xmax(0), int(self._Nt / 2)) - stopnumb + 1 ]
        else:
            self._fit_interval[1] += 1




    def xmax(self, i = None):
        if i is None:
            if self._fit_interval is None:
                return self.xmax(self._nt_min_start)
            else:
                return self.xmax(self._fit_interval[0])
        if self._xmax is None:
            if self._sym:
                if self._ratio_fit:
                    return int(self._Nt / 2 - math.floor(self._Nt / 12))
                else:
                    return int(self._Nt / 2)
            else:
                return int(self._Nt - i)
        else:
            return self._xmax


    def xmin(self):
        if self._fit_interval is not None:
            return self._fit_interval[0]
        else:
            return self._nt_min_start


    def _get_start_params(self, res, res_err, ratio = np.inf, ignore_ext = False):

        if not ignore_ext and self._start_params is not None:
            logger.info("Use start parameters from user")
            return self._start_params

        if len(res) < 1 or len(res_err) < 1:
            return None

        last_res = res[-1]
        last_err = res_err[-1]

        for i in range(len(last_res)):

            if np.isnan(last_res[i]):
                return None

            if last_err[i] / np.abs(last_res[i]) > ratio:
                return None

        return last_res



    def check_data_consistency(self):

        if len(self._xdata) != len(self._data):
            raise ValueError("Length of xdata != lengths of data")

        if self._Nt%2 != 0:
            raise ValueError("Nt must be even!")

        last_len = len(self._data[0])
        for i in self._data:
            if len(i) != last_len:
                raise ValueError("Number of data points per x value inconsistent!")


    def read_in_pure(self):
        symmetrize = (self._Nt is None) & self._sym
        self._xdata, self._data, self._nconfs = read_in_pure(self._filename, self._nt_col,
                self._data_col, symmetrize = symmetrize)

        if self._Nt is None:
            self._Nt = len(self._data)

        self.check_data_consistency()

        self._calc_fit_interval()

        ind = (self._xdata <= self.xmax()) & (self._xdata >= self.xmin())
        self._xdata = self._xdata[ind]
        self._data = self._data[ind]

        if self._auto_sign:
            mean = np.mean(self._data, axis = 1)
            even = self._xdata % 2 == 0
            test_neg = np.sum(0 > mean[even])
            test_pos = np.sum(0 < mean[even])
            if test_neg > test_pos:
                self._change_sign = True

        if self._change_sign:
            logger.warn("\n\nCHANGING SIGN OF CORRELATOR.\nThis interchanges oscillating and "
                    "non-oscillating part!\n")
            for i in range(len(self._data)):
                self._data[i] *= (-1)**(self._xdata[i]+1)






    def read_in_jk(self):
        self.read_in_pure()

        self._ydata, self._edata = jackknife(
                self._calc_mean, self._data, self._numb_blocks)

        if self._correlated:
            self._cov = self._calc_cov()




    def read_in_std_err(self):
        self.read_in_pure()

        if self._correlated:
            self._ydata, self._edata, self._cov = self._calc_cov_and_mean()
        else:
            self._ydata, self._edata = mean_and_err(self._data, axis = 1)



    def read_in_sample(self):
        self.read_in_pure()

        if self._correlated:
            self._ydata, self._edata, self._cov = self._calc_cov_and_mean()
        else:
            self._ydata, self._edata = mean_and_std_dev(self._data, axis = 1)



    def read_in_direct(self):
        self._xdata, self._ydata, self._edata = read_in(self._filename,
                self._nt_col, self._data_col, self._err_col)
        if self._Nt is None:
            self._Nt = len(self._ydata)
        self._calc_fit_interval()



    def read_in_jk_ratio(self):
        self.read_in_jk()
        self._ratio_ydata, self._ratio_edata = jackknife(corr_ratio,
                self._data, self._numb_blocks)
        self._ratio_xdata = self._xdata[:-1]
        if self._correlated:
            self._ratio_cov = self._calc_ratio_cov()




    def read_in_std_err_ratio(self):
        self.read_in_std_err()
        self._rat_data = np.array([ corr_ratio_direct(i)
            for i in self._data.transpose()]).transpose()

        self._ratio_ydata, self._ratio_edata = mean_and_err(self._rat_data, axis=1)
        self._ratio_xdata = self._xdata[:-1]
        if self._correlated:
            self._ratio_cov = self._calc_ratio_cov()



    def read_in_bs(self):
        self.read_in_pure()
        # If we use the correlation matrix, that ydata should keep their correlation. Therefore
        # same_rand_for_obs = self._correlated

        self._ydata, self._edata = bs.bootstr(
                self._calc_mean, self._data, self._numb_samples,
                same_rand_for_obs = self._correlated,
                err_by_dist = self._ng_btstr, seed = self._seed)

        if self._correlated:
            self._cov = self._calc_cov()





    def read_in_results(self):
        (self._ranges, self._res, self._res_err, self._aicc,
                self._chi_dof, self._nstates, self._nstates_osc) = read_in_fitmass(
                self._res_filename)

        self._nparams = 2 * (self._nstates + self._nstates_osc)



    def _calc_mean(self, data = None):

        if data is None:
            data = self._data

        if self._min_cov_det:
            logger.details("Use MCD for expectation value estimation")
            mcd = MCD(support_fraction = self._mcd_supp_frac).fit(self._data.transpose())
            ydata = mcd.location_
        else:
            ydata = np.mean(data, axis = 1)

        return ydata




    def _calc_cov(self, data = None):

        if data is None:
            data = self._data
        if self._min_cov_det:
            logger.details("Use MCD for covariance estimation")
            mcd = MCD(support_fraction = self._mcd_supp_frac).fit(self._data.transpose())
            cov = mcd.covariance_
        else:
            cov = calc_cov(data.transpose())

        if not self._sample_data:
            cov /= self._nconfs  # For fit we have to normalize like an error
        return cov



    def _calc_cov_and_mean(self, data = None):

        if data is None:
            data = self._data
        if self._min_cov_det:
            ind = self._xdata < self.xmax()
            mcd = MCD(support_fraction = self._mcd_supp_frac).fit(self._data.transpose())
            ydata = mcd.location_
            logger.details("Use MCD for covariance estimation")
            cov = mcd.covariance_
        else:
            ydata = std_mean(self._data, axis = 1)
            cov = calc_cov(data.transpose())

        if not self._sample_data:
            cov /= self._nconfs  # For fit we have to normalize like an error
        edata = np.sqrt(np.diag(cov))
        return ydata, edata, cov


    def _calc_ratio_cov(self, data = None):
        return self._calc_cov(self._rat_data)





    """Compute the effective mass for all available distances.
    In case of nstates_ocs > 0 we
    get two different effective masses for the oscillating and non-oscillating part.
    (See https://arxiv.org/abs/1411.3018).
    In that case we also compute the split correlators. These are needed for the parameter
    estimation for an oscillating fit."""
    def _compute_meff(self):
        logger.info("Cmpute effective mass...")


        # It is always better to use a jackknife for meff instead of a computation
        # based on error propagation. Therefore we use a Jacknife also if the data flag is set
        # to std_err_data.
        if self._jack_data or self._std_err_data or self._ratio_data:
            if self._nstates_osc > 0:
                rm_osc, rm_osc_err = jackknife(
                        remove_osc, self._data, self._numb_blocks, args=(self._xdata,))

            else:
                eff_mass, eff_mass_err = jackknife(calc_eff_mass, self._data,
                        self._numb_blocks, args = (self._xdata, self._Nt))



        if self._btstr_data:

            if self._nstates_osc > 0:
                rm_osc, rm_osc_err = bs.bootstr(remove_osc, self._data, self._numb_samples,
                        err_by_dist = self._ng_btstr, args=(self._xdata,),
                        seed = self._seed)
            else:
                eff_mass, eff_mass_err = bs.bootstr(calc_eff_mass, self._data,
                        self._numb_samples, err_by_dist = self._ng_btstr,
                        args = (self._xdata, self._Nt), seed = self._seed)



        if self._sample_data:
            if self._nstates_osc > 0:
                results = []
                for i in range(len(self._data[0])):
                    results.append(remove_osc_av(self._data[:, i], self._xdata))
                rm_osc = std_mean(results)
                rm_osc_err = std_dev(results)

            else:
                eff_masses = []
                for i in range(len(self._data[0])):
                    eff_masses.append(calc_eff_mass_direct(self._data[:, i], self._xdata,
                        self._Nt))
                eff_mass, eff_mass_err = mean_and_std_dev(eff_masses)



        if not self._direct_data:
            if self._nstates_osc > 0:

                # Nan destroys the whole fit, as the chi^2 becomes nan. Therefore we have to remove
                # those from the helper arrays for parameter estimation.
                rm_osc = remove_nan(*(rm_osc + rm_osc_err))

                self._xdata_osc = rm_osc[0]
                self._xdata_no = rm_osc[0]
                self._ydata_no = rm_osc[1]
                self._ydata_osc = rm_osc[2]
                self._edata_no = rm_osc[8]
                self._edata_osc = rm_osc[9]


                self._meff_no = rm_osc[3]
                self._meff_osc = rm_osc[4]
                self._meff_no_err = rm_osc[10]
                self._meff_osc_err = rm_osc[11]


            else:
                self._meff = eff_mass
                self._meff_err = eff_mass_err





    def init_fitter(self):

        if self._correlated:
            edata = self._cov
        else:
            edata = self._edata

        if self._ratio_fit:
            if self._correlated:
                self._corr_fitter = CorrFitterRatio(self._ratio_xdata, self._ratio_ydata,
                        self._ratio_cov, self._xdata, self._ydata, self._edata,
                        self._nstates, Nt = self._Nt)
            else:
                self._corr_fitter = CorrFitterRatio(self._ratio_xdata, self._ratio_ydata,
                        self._ratio_edata, self._xdata, self._ydata, self._edata,
                        self._nstates, Nt = self._Nt)
        else:
            self._corr_fitter = CorrFitter(self._xdata, self._ydata, edata,
                    nstates = self._nstates, nstates_osc = self._nstates_osc,
                    cut_eig = self._cut_eig, cut_perc = self._cut_perc, Nt = self._Nt)

        if self._nstates_osc > 0:
            # Helper fitter for parameter estimation
            if len(self._ydata_osc) > 0:
                self._osc_corr_fitter = CorrFitter(self._xdata_osc, self._ydata_osc,
                        self._edata_osc, Nt = self._Nt,
                        nstates = self._nstates_osc, nstates_osc = 0)

                # We find that for the split correlators, excited states do not contribute
                # that much in the middle of the correlator. Therefore we set this weight to
                # a higher value
                self._osc_corr_fitter.est_weights[2] = 4

            if len(self._ydata_no) > 0:
                self._no_corr_fitter = CorrFitter(self._xdata_no, self._ydata_no,
                        self._edata_no, Nt = self._Nt,
                        nstates = self._nstates, nstates_osc = 0)

                # We find that for the split correlators, excited states do not contribute
                # that much in the middle of the correlator. Therefore we set this weight to
                # a higher value
                self._no_corr_fitter.est_weights[2] = 4




    def try_fit(fit_func):
        def func_wrapper(inst, *args, **kwargs):
            try:
                res, res_err, chi_dof, aicc, pcov = fit_func(inst, *args, **kwargs)
                if any(res_err == 0.0):
                    raise ValueError("Fit error is zero!")
            except Exception as e:
                logger.warn("\nFit failed! Exception was:", e, "\n")
                if logger.isLevel("DEBUG"):
                    traceback.print_exc()
                res = np.full(inst._nparams, np.nan)
                res_err = np.full(inst._nparams, np.nan)
                pcov = np.full((inst._nparams, inst._nparams), np.nan)
                chi_dof = np.inf
                aicc = np.inf
                if inst._btstr_fit:
                    inst._samples.append([[(np.nan,)*inst._nparams, (np.nan,)*inst._nparams,
                        np.nan, np.nan] for sitt in range(inst._nfit_samples)])
            return res, res_err, chi_dof, aicc, pcov
        return func_wrapper





    @try_fit
    def direct_fit(self, xmin, xmax, start_params):
        if self._ratio_fit:
            return self._corr_fitter.corr_fit(xmin, xmax, start_params,
                    priorval = self._priorval, priorsigma = self._priorsigma)
        else:
            return self._corr_fitter.corr_fit(xmin, xmax, start_params,
                    nstates = self._nstates, nstates_osc = self._nstates_osc,
                    priorval = self._priorval, priorsigma = self._priorsigma,
                    correlated = self._correlated)





    def _fit_corr_fr_data(self, data, xmin, xmax, start_params):

        logger.info()
        logger.info("Fitrange [ %d, %d ]" %(xmin, xmax))

        if self._scnd_btstr:
            ydata, edata = bs.bootstr(
                    self._calc_mean, data, self._numb_samples,
                    err_by_dist = self._ng_btstr, same_rand_for_obs = False)
            corr_fitter = CorrFitter(self._xdata, ydata, edata, Nt = self._Nt)
        else:
            if self._correlated:
                ydata, edata, cov = self._calc_cov_and_mean(data)
                corr_fitter = CorrFitter(self._xdata, ydata, cov, Nt = self._Nt)

            else:
                ydata, edata = mean_and_err(data, axis=1)
                corr_fitter = CorrFitter(self._xdata, ydata, edata, Nt = self._Nt)

        res, res_err, chi_dof, aicc, pcov = corr_fitter.corr_fit(xmin, xmax, start_params,
                correlated = self._correlated, priorval = self._priorval,
                priorsigma = self._priorsigma, nstates = self._nstates,
                nstates_osc = self._nstates_osc)

        if self._nstates_osc > 0:
            if res[1] < res_err[1] or res[2*self._nstates + 1] < res_err[2*self._nstates + 1]:
                raise ValueError("Fit error larger than value")
        else:
            if res[1] < res_err[1]:
                raise ValueError("Fit error larger than value")

        return res, res_err, chi_dof, aicc, pcov




    @try_fit
    def jk_fit(self, xmin, xmax, start_params):
        logger.info("\nFirst fit to estimate start parameters...")
        start_params = self.direct_fit(xmin, xmax, start_params)[0]
        logger.info("\nStart jackknife...")

        ((res, tmp_err, chi_dof, aicc, pcov),
                (res_err, tmp_err_err, chi_err, aicc_err, pcov_err)) = jackknife(
                    self._fit_corr_fr_data, self._data, self._numb_blocks,
                    args = (xmin, xmax, start_params) )

        logger.info("\n\n\n")
        print_res("Final result for %d + %d jackknife fit"
                % (self._nstates, self._nstates_osc), res, res_err, chi_dof, level = "INFO")

        return res, res_err, chi_dof, aicc, pcov


    @try_fit
    def bs_fit(self, xmin, xmax, start_params):
        logger.info("\nFirst fit to estimate start parameters")

        start_params = self.direct_fit(xmin, xmax, start_params)[0]
        logger.info("\nStart bootstrap...")

        if self._correlated:
            (sample, (res, tmp_err, chi_dof, aicc, pcov),
                    (res_err, tmp_err_err, chi_err, aicc_err, pcov_err)) = bs.bootstr(
                    self._fit_corr_fr_data, self._data, self._nfit_samples,
                    args = (xmin, xmax, start_params), same_rand_for_obs = True,
                    err_by_dist = self._ng_btstr,
                    return_sample = True, seed = self._seed,
                    nmax_exceptions = 0.5 * self._nfit_samples)
        else:
            (sample, (res, tmp_err, chi_dof, aicc, pcov),
                    (res_err, tmp_err_err, chi_err, aicc_err, pcov_err)) = bs.bootstr(
                    self._fit_corr_fr_data, self._data, self._nfit_samples,
                    args = (xmin, xmax, start_params), return_sample = True,
                    err_by_dist = self._ng_btstr, seed = self._seed,
                    nmax_exceptions = 0.5 * self._nfit_samples)


        self._samples.append(sample)
        all_res = [ i[0] for i in sample ]

        # Compute covariance matrix from sample
        pcov = calc_cov(np.array(all_res).transpose())

        logger.info("\n\n\n")
        print_res("Final result for %d + %d bootstrap fit"
                % (self._nstates, self._nstates_osc), res, res_err, chi_dof, level = "INFO")


        return res, res_err, chi_dof, aicc, pcov





    # Get the starting parameters for a non-oscillating fit.
    # This is only possible if the original data is available
    def _get_start_params_osc(self, xmin, xmax):
        self._osc_fit(xmin, xmax)

        self._no_fit(xmin, xmax)

        self._add_sep_data()

        if any(np.isnan(self._sep[-1])):
            start_params = None
            logger.info("Failed to estimate start parameters for oscillating fit")
        else:
            start_params = self._sep[-1]

        return start_params


  # oscillating part
    def _osc_fit(self, xmin, xmax):
        # We do not use constraint fits here as this fit is for parameter estimation
        logger.info()
        start_params_osc = self._get_start_params(self._res_osc, self._res_osc_err,
                ratio = 1, ignore_ext = True)

        try:
            logger.info("Fit osc. correlator ...")
            (res_osc, res_osc_err, chi_dof_osc, aicc_osc,
                    pcov_osc) = self._osc_corr_fitter.corr_fit(
                        xmin, xmax, start_params_osc, correlated = False,
                        nstates = self._nstates_osc, nstates_osc = 0)

            self._res_osc.append(res_osc)
            self._res_osc_err.append(res_osc_err)
            self._chi_dof_osc.append(chi_dof_osc)
            self._aicc_osc.append(aicc_osc)

        except Exception as e:
            logger.warn(e, "\n")
            if logger.isLevel("DEBUG"):
                traceback.print_exc()

            self._res_osc.append(np.full(2 * self._nstates_osc, np.nan))
            self._res_osc_err.append(np.full(2 * self._nstates_osc, np.nan))
            self._chi_dof_osc.append(np.nan)
            self._aicc_osc.append(np.nan)




    # non-oscillating part
    def _no_fit(self, xmin, xmax):
        # We do not use constraint fits here as this fit is for parameter estimation
        logger.info()
        start_params_no = self._get_start_params(self._res_no, self._res_no_err,
                ratio = 1, ignore_ext = True)

        try:
            logger.info("Fit non-osc. correlator...")
            res_no, res_no_err, chi_dof_no, aicc_no, pcov_no = self._no_corr_fitter.corr_fit(
                    xmin, xmax, start_params_no, correlated = False,
                    nstates = self._nstates, nstates_osc = 0)


            self._res_no.append(res_no)
            self._res_no_err.append(res_no_err)
            self._chi_dof_no.append(chi_dof_no)
            self._aicc_no.append(aicc_no)

        except Exception as e:
            logger.warn("\nFit failed:", e, "\n")
            if logger.isLevel("DEBUG"):
                traceback.print_exc()

            self._res_no.append(np.full(2 * self._nstates, np.nan))
            self._res_no_err.append(np.full(2 * self._nstates, np.nan))

            self._chi_dof_no.append(np.nan)
            self._aicc_no.append(np.nan)



    # Store data of osc. and non-osc. fit as if they come from a combined fit
    def _add_sep_data(self):
        sep = list(self._res_no[-1]) + list(self._res_osc[-1])
        sep_err = list(self._res_no_err[-1]) + list(self._res_osc_err[-1])
        self._sep.append(sep)
        self._sep_err.append(sep_err)

        # This does not make sense, but we have to fill the arrays somehow.
        self._chi_dof_sep.append(self._chi_dof_no[-1] + self._chi_dof_osc[-1])
        self._aicc_sep.append(self._aicc_no[-1] + self._aicc_osc[-1])





    def _get_start_params_array(self, xmin, xmax):

        start_params_array = []

        if self._start_params is None:


            # If we try all methods, we want to use the parameters of the last fit.
            # No matter how bad they are. Therefore set ratio to inf.
            if self._try_all:
                ratio = np.inf
            else:
                ratio = 1

            start_params = self._get_start_params(self._res, self._res_err, ratio = ratio)

            if start_params is not None:
                start_params_array.append((start_params, "last fit"))

            if not (self._osc_corr_fitter is None
                    or self._no_corr_fitter is None) and not self._ratio_fit:

                start_params = self._get_start_params_osc(xmin, xmax)

                if self._try_all or len(start_params_array) == 0:
                    if start_params is not None:
                        start_params_array.append((start_params, "separate fits"))


            if self._try_all or len(start_params_array) == 0:
                start_params_array.append((None, "data"))
        else:
            start_params_array.append((self._start_params, "user"))

        return start_params_array


    def perform_fit(self, xmin, xmax):
            logger.info("==========================================================================")
            logger.info("Fitrange [ %d, %d ]" %(xmin, xmax))

            self._ranges.append([xmin, xmax])


            start_params_array = self._get_start_params_array(xmin, xmax)

            all_res = []
            all_chi_dof = []
            for (start_params, text) in start_params_array:
                logger.info("\n--------------------------------------")
                logger.info("Run with start parameters from", text, "\n")
                res = self.fit(xmin, xmax, start_params)
                all_chi_dof.append(res[2])
                all_res.append(res)

            # Find smallest chi^2
            min_ind = np.argmin(all_chi_dof)
            logger.info("\n--------------------------------------")
            logger.info("Choose fit with start parameters from",
                    start_params_array[min_ind][1], "\n")


            res, res_err, chi_dof, aicc, pcov = all_res[min_ind]

            print_res("Final fit_res for xmin = %d, xmax = %d" % (xmin, xmax),
                    res, res_err, chi_dof, level = "INFO")
            print_scl("AICc", aicc, level = "INFO")

            self._res.append(res)
            self._res_err.append(res_err)
            self._pcov.append(pcov)
            self._chi_dof.append(chi_dof)
            self._aicc.append(aicc)

            logger.info()
            return res, res_err, chi_dof, aicc



    def _fit_loop(self):
        logger.info("\nStart fit loop from tmin = %d to tmax = %d"
                % (self._fit_interval[0], self._fit_interval[1] - 1))
        logger.info("\nYou will probably see a lot of errors\n"
                "This is normal, as we try a lot of different fitting methods\n"
                "and not all of them will work\n")
        for i in range(*self._fit_interval):
            xmin, xmax = i, self.xmax(i)
            self.perform_fit(xmin, xmax)






    def _plot_corr_wrapper(self, fitter, res, name, add_to_title, no_error = False):
        if self._title is not None:
            title = add_to_title + self._title
        else:
            title = None

        os.makedirs(self._folder, exist_ok=True)
        fitter.plot_corr(
                self._folder + "/corr_" + name + self._out_name + ".pdf",
                res, None, self._ranges,
                notex = self._notex, title = title,
                xlabel = self._xlabel, ylabel = self._ylabel, size = self._plot_size,
                font_size = self._font_size, plot_ylog = self._plot_ylog,
                no_error = no_error
                )



    def plot_corr(self):
        # Use the covariance matrix of the parameters for computing the error bands via error
        # propagation
        fitter = self._corr_fitter
        self._plot_corr_wrapper(fitter, self._res, "", "",
                no_error = True)

        if self._nstates_osc > 0:

            if len(self._sep) > 0:
                self._plot_corr_wrapper(self._corr_fitter, self._sep,
                        "sep_", "sep., ", no_error = True)

            if len(self._res_osc) > 0:
                self._plot_corr_wrapper(self._osc_corr_fitter, self._res_osc,
                        "osc_", "osc., ", no_error = True)

            if len(self._res_no) > 0:
                self._plot_corr_wrapper(self._no_corr_fitter, self._res_no,
                        "no_", "non-osc., ", no_error = True)




    def _plot_cov_mat(self):
        os.makedirs(self._folder, exist_ok=True)


        xmin = np.min(self._fit_interval)
        xmax = self.xmax(xmin)

        self._corr_fitter.plot_cov(filename = self._folder + "/cov_" + self._out_name + ".pdf",
                title = self._title, notex = self._notex,
                xmin = xmin, xmax = xmax,
                ymin = xmin, ymax = xmax)

        self._corr_fitter.plot_eig(filename = self._folder + "/eig_" + self._out_name + ".pdf",
                title = self._title, notex = self._notex, xmin = xmin, xmax = xmax)


    def _write_data(self):
        os.makedirs(self._folder, exist_ok=True)

        if self._sym:
            method_string = "sym"
        else:
            method_string = "asym"

        if self._jack_data or self._std_err_data or self._ratio_data:
            method_string += "_jk"

        if self._btstr_data:
            method_string += "_btstr"

        if self._sample_data:
            method_string += "_fr_sample"


        if not self._direct_data:
            logger.info("Write effective mass...")
            if self._nstates_osc > 0:
                write_eff_mass(self._folder + "/effmass_no_" + method_string
                        + self._file_string + ".txt", self._meff_no, self._meff_no_err,
                        self._xdata_no, osc = True)

                write_eff_mass(self._folder + "/effmass_osc_" + method_string
                        + self._file_string + ".txt", self._meff_osc, self._meff_osc_err,
                        self._xdata_osc, osc = True)
            else:
                write_eff_mass(self._folder + "/effmass_" + method_string
                        + self._file_string + ".txt", self._meff, self._meff_err, self._xdata)


        write_fit_mass(self._folder + "/fitmass_" + self._out_name + ".txt",
                self._ranges, self._res, self._res_err, self._chi_dof, self._aicc,
                self._nstates, self._nstates_osc)
        write_pcov(self._folder + "/pcov_" + self._out_name + ".txt", self._ranges, self._pcov,
                self._nstates, self._nstates_osc)
        write_corr(self._folder + "/corr_" + self._out_name + ".txt",
                self._xdata, self._ydata, self._edata)
        if self._btstr_fit:
            write_sample(self._folder + "/fitmass_sample_" + self._out_name + ".txt",
                    self._ranges, self._samples, self._nstates, self._nstates_osc)


    def get_results(self):
        return self._res, self._res_err, self._chi_dof, self._aicc

    def get_xyedata(self):
        return self._xdata, self._ydata, self._edata

    def get_data(self):
        return self._data

    def get_cov(self):
        return self._cov

    def get_Nt(self):
        return self._Nt

    def get_nstat(self):
        return len(self._data[0])

    def run(self):

        if self._res_filename is not None:
            self._out_name += "_fr-file"
            self.read_in_results()



        if not self._plot_start:
            if self._res_filename is None:
                self._fit_loop()

        else:
            self._res.append(self._start_params)
            fit_range = (self._fit_interval[0], self.xmax(self._fit_interval[0]))
            self._ranges.append(fit_range)
            self._res_err.append(np.zeros_like(self._start_params))
            self._out_name += "_start_params"
            print_scl("Chi^2/d.o.f.",
                    self._corr_fitter.calc_chisquare_dof(self._res[-1], *fit_range),
                    level = "INFO")



        logger.info("Plot data...")
        self.plot_corr()
        if self._correlated:
            logger.info("Plot covariance matrix...")
            self._plot_cov_mat()
        if not (self._plot_start or self._res_filename is not None):
            logger.info("Write out data...")
            self._write_data()
