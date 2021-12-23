import numpy as np
from latqcdtools.fitting import Fitter, print_res, print_scl, NotAvailableError
from latqcdtools.tools import numToWords, merge_two_dicts
from scipy.optimize import curve_fit
import latqcdtools.logger as lg
import traceback


class CorrFitter(Fitter):

    # If edata contains errors, they are not squared: edata=sqrt(variance). If it is the covariance matrix they are
    # interpreted as squared:
    def __init__(self, xdata, ydata, edata = None,
            Nt = None, nstates = 1, nstates_osc = 0, **kwargs):
        if Nt is None:
            Nt = len(xdata)
        self._Nt = Nt

        self.est_weights = [1, 1, 2, 4, 6, 6, 6, 12, 14]
        self.est_weights_osc = [1, 1, 8, 8, 10, 10, 10, 12, 14]
        # To make sure we have floats.
        xdata = np.array(xdata, dtype = float)

        # For convience we also save error data
        if edata is not None:
            try:
                edata[0][0]
                self._edata = np.sqrt(np.diag(edata))
            except (IndexError, TypeError) as e:
                self._edata = np.asarray(edata)
        else:
            self._edata = np.ones(len(ydata))

        # Initialize variables for the fitting function before the initialization of the fitter
        self._nstates = nstates
        self._nstates_osc = nstates_osc
        self._mult_const = 1

        if kwargs.get('cut_eig', False):
            use_corr = False
        else:
            use_corr = True


        Fitter.__init__(self, self._corr, xdata, ydata, edata,
                grad = self._grad_corr, hess =  self._hess_corr, expand = False,
                args = (self._nstates, self._nstates_osc), use_corr = use_corr,
                func_sup_numpy = True, **kwargs)

        self.calc_mult_const()

        # This defines the method that is used for the fit. It is still possible, to fit with
        # different methods via passing a different enum to simple_corr_fit
        self._fit_name = ""


  

    # It turns out that the fit works much better, if we extract a multiplication constant,
    # such that the amplitude is roughly 1.0
    def calc_mult_const(self):
        try:
            mult_const = np.abs(self._ydata[int(self._Nt / 2)])
        except IndexError:
            try:
                mult_const = abs(self._ydata[-1])
            #If data are empty
            except IndexError:
                mult_const = 1


        lg.info("Normalize correlator with", mult_const)
        self._mult_const = abs(mult_const)
        return self._mult_const



    def set_mult_const(self, mult_const):
        self._mult_const = mult_const

    def get_mult_const(self):
        return self._mult_const


    def set_states(self, nstates, nstates_osc):
        self._nstates = nstates
        self._nstates_osc = nstates_osc
        self._args = (nstates, nstates_osc)
        self._grad_args = (nstates, nstates_osc)
        self._hess_args = (nstates, nstates_osc)


    """
    For higher state fits, it might happen, that the order of the results does not correspond
    to the order of the states. Therefore we have to reorder the results
    This is done by looping over the results and checking, whether two states need to be
    exchanged.
    We exchange the two states, if
    - the amplitude of the first state is lower than the second one's
    - the mass of the first state is unrealistic small ( < 1e-5)
    - the amplitude of the first state is negative and the amplitude of the second state is
      positive
    """

    def _swap(self, i, res, res_err):
        if ((abs(res[i + 2]) > abs(res[i + 0]) and abs(res[i + 3]) > 1e-5
                or (abs(res[i + 1]) < 1e-5 < abs(res[i + 3]))
                or (res[i + 0] < 0 < res[i + 2]))
                and not (res[i + 0] > 0 > res[i + 2])):

            res[i + 1], res[i + 3] = abs(res[i + 3]), abs(res[i + 1]),
            res[i + 0], res[i + 2] = res[i + 2], res[i + 0]

            res_err[i + 1], res_err[i + 3] = res_err[i + 3], res_err[i + 1],
            res_err[i + 0], res_err[i + 2] = res_err[i + 2], res_err[i + 0]

    """
    This is the loop to change the order of the states. 
    Loops for oscillating and non-oscillating case are seperated
    """
    def _change_order(self, res, res_err, nstates, nstates_osc): 

        for j in range(0, nstates):
            for i in range(0, 2*nstates - 3, 2):
                self._swap(i, res, res_err)

        for j in range(0, nstates_osc):
            # Second state data start at 2*nstates
            for i in range(2 * nstates, 2 * (nstates + nstates_osc) - 3, 2):
                self._swap(i, res, res_err)


    """A simple correlator fit without starting parameter estimation"""
    def simple_corr_fit(self,
            xmin = -np.inf, xmax = np.inf,
            start_params = None, correlated = None,
            algorithms = None,
            nstates = 1, nstates_osc = 0, 
            priorsigma = None, priorval = None,
            ):

        if algorithms is None:
            algorithms = CorrFitter._std_algs

        self._args = (nstates, nstates_osc)
        self._grad_args = (nstates, nstates_osc)
        self._hess_args = (nstates, nstates_osc)

        res, res_err, chi_dof, pcov = self.try_fit(start_params = start_params, xmin = xmin, xmax = xmax,
                algorithms = algorithms, correlated = correlated,
                priorsigma = priorsigma, priorval = priorval, ret_pcov = True)
        aicc = self.aicc(res, xmin, xmax, correlated = correlated)
        return res, res_err, chi_dof, aicc, pcov


    def remove_mult_const(self, res, res_err):
        for i in range(0,len(res), 2):
            res[i] *= self._mult_const
            res[i+1] = abs(res[i+1])
            res_err[i] *= self._mult_const
            res_err[i+1] = abs(res_err[i+1])
        return res, res_err


    def _apply_mult_const(self, res):
        for i in range(0,len(res), 2):
            res[i] /= self._mult_const
        return res



    def remove_mult_const_pcov(self, pcov):
        for j in range(0, len(pcov)):
            for i in range(0, len(pcov[j])):
                if i % 2 == 0:
                    pcov[i,j] *= self._mult_const
                if j % 2 == 0:
                    pcov[i,j] *= self._mult_const
        return pcov
 
 
    def init_start_params(self,
            xmin, xmax, start_params, priorval, priorsigma, nstates, nstates_osc
            ):

        if priorval is not None:
            start_params = np.copy(priorval)
            priorval = np.copy(priorval)
            priorsigma = np.copy(priorsigma)
            for i in range(0,len(priorval),2):
                priorval[i] /= self._mult_const
                priorsigma[i] *= self._mult_const

        if start_params is not None:
            start_params = np.copy(start_params)
            for i in range(0,len(start_params),2):
                start_params[i] /= self._mult_const
        else:
            start_params = self._est_params(xmin, xmax, nstates, nstates_osc)

        return start_params, priorval, priorsigma


  
    def corr_fit(self,
            xmin = -np.inf, xmax = np.inf, start_params = None, 
            nstates = None, nstates_osc = None, correlated = None,
            priorsigma = None, priorval = None):


        if nstates == 0:
            raise ValueError("Require at least one non-oscillating state")

        if nstates is None:
            nstates = self._nstates
        if nstates_osc is None:
            nstates_osc = self._nstates_osc

        if correlated is None:
            correlated = self._cov_avail

        # To estimate parameters we usually perform a non-correlated fit before the actual
        # correlated fit. This is not neccessary if we already get start parameters.
        # However, for an oscillating fit is is reasonable to perform a non correlated fit
        # in any case. This is because start parameter estimation is usually performed outside
        # for oscillating fits.
        if correlated and start_params is not None:
            if nstates_osc > 0:
                skip_uncorr = False
            else:
                skip_uncorr = True
        else:
            skip_uncorr = False

        try:
            start_params, priorval, priorsigma = self.init_start_params(
                    xmin, xmax, start_params, priorval, priorsigma, nstates, nstates_osc)
        except Exception as e:
            lg.info("Failed to estimate start parameters. Try direct fit")
            lg.details("Error was", e)
            if lg.isLevel("DEBUG"):
                traceback.print_exc()
            start_params = None

        # Save the states of the last fit. This must be ensured, even if the parameter
        # estimation fails.
        finally:
            self._nstates = nstates
            self._nstates_osc = nstates_osc


        print_res("Start parameters for %d + %d fit" % (nstates, nstates_osc),
                start_params, level = "INFO")

        if not skip_uncorr:
            res, res_err, chi_dof, aicc, pcov = self.simple_corr_fit(xmin, xmax, start_params,
                    correlated = False, priorval = priorval, priorsigma = priorsigma, 
                    nstates = nstates, nstates_osc = nstates_osc
                    )

            start_params = np.copy(res)


            self._change_order(res, res_err, nstates, nstates_osc)
            res, res_err = self.remove_mult_const(res, res_err)
            pcov = self.remove_mult_const_pcov(pcov)

            lg.info()
            print_res("Fit result for uncorrelated %d + %d fit" % (nstates, nstates_osc),
                    res, res_err, chi_dof, level = "INFO")

            print_scl("AICc", aicc, level = 'INFO')

        if correlated:
            if not self._cov_avail:
                raise NotAvailableError("Covariance matrix is not available")

            res, res_err, chi_dof, aicc, pcov = self.simple_corr_fit(
                    xmin, xmax, start_params = start_params, correlated = True,
                    priorval = priorval, priorsigma = priorsigma,
                    nstates = nstates, nstates_osc = nstates_osc
                    )


            self._change_order(res, res_err, nstates, nstates_osc)
            res, res_err = self.remove_mult_const(res, res_err)
            pcov = self.remove_mult_const_pcov(pcov)

            lg.info()
            print_res("Fit result for correlated %d + %d fit" % (nstates, nstates_osc),
                    res, res_err, chi_dof, level = "INFO")

            print_scl("AICc", aicc)

        return res, res_err, chi_dof, aicc, pcov



  

    def _get_est_range(self, xmin, xmax, weight_low = 2, weight_up = 1):
        if xmax > int(self._Nt / 2):
            est_range = [int((weight_low * xmin + weight_up * self._Nt / 2)
                / (weight_up + weight_low)),
                int((weight_low * (xmax + 1) + weight_up * self._Nt / 2) /
                    (weight_up + weight_low))]
        else:
            est_range = [int((weight_low * xmin + weight_up * xmax) /
                (weight_up + weight_low)), xmax]
        return est_range


    def _est_params(self, xmin, xmax, nstates, nstates_osc):
        try:
            if nstates_osc > 0:
                return self._est_params_osc(xmin, xmax, nstates, nstates_osc)
            else:
                return self._est_params_non_osc(xmin, xmax, nstates)
        finally:
            #In case the above fit failed, we have to reset the function
            self._func = self._corr
            self._args = (nstates, nstates_osc)
            self._grad_args = (nstates, nstates_osc)
            self._hess_args = (nstates, nstates_osc)


    """Estimate start parameters for one state fit without oscillations"""
    def _est_params_one_state(self, xmin, xmax, ind = None):

        # estimation of starting parameters via linear fit does not work for periodic
        # boundaries. Therefore we have to fit on the half interval.
        start_params = [1, 1]
        est_range = [int(xmin), int(min(xmax, self._Nt/2))]

        # For estimating parameters of oscillating correlators, we use this function to estimate
        # the parameters for the even/odd part. We have to make sure that we use even/odd part
        # by this routine.
        lg.info("Estimate parameters: Linear fit for one state...")
        if ind is None:
            ind = (self._xdata > est_range[0]) & (self._xdata < est_range[1])
        else:
            ind = (self._xdata > est_range[0]) & (self._xdata < est_range[1]) & ind

        try:
            (tmp_start_params, tmp_start_params_cov) = curve_fit(linear, self._xdata[ind],
                    np.log(np.abs(self._ydata[ind])), sigma=1 / self._ydata[ind] * self._edata[ind])

            if tmp_start_params[0] > 50 or tmp_start_params[1] > 50:
                start_params = [1, 1]
            else:
                start_params[0] = 1 / self._mult_const * np.exp(tmp_start_params[0])\
                    * np.exp(-tmp_start_params[1] * self._Nt / 2) * 2

                start_params[1] = tmp_start_params[1]
        except Exception as e:
            lg.info("Linear fit failed. Error was", e)
            if lg.isLevel("DEBUG"):
                traceback.print_exc()
            start_params = [1, 1]

        print_res("Final start parameters for uncorrelated fit", start_params, level = "INFO")

        return start_params



    """Estimate start parameters for multiple state fit without oscillations"""
    def _est_params_non_osc(self, xmin, xmax, nstates):
        lg.info("Estimate parameters: Fit range for", numToWords(nstates),
                "state fit: [", xmin, ",", xmax, "]")
        nparams = 2 * nstates
        start_params = np.ones(nparams)

        est_range = self._get_est_range(xmin, xmax, self.est_weights[nstates])


        if nstates == 1:
            return self._est_params_one_state(xmin, xmax)
        else:
            start_params[0:nparams - 2] = self._est_params_non_osc(est_range[0], est_range[1],
                    nstates - 1)

        print_res("Starting parameters for " + numToWords(nstates - 1) + " state fit",
                start_params[0:nparams - 2], level = "INFO")

        start_params[0:nparams - 2] = self.simple_corr_fit(
                nstates = nstates - 1, correlated = False, xmin = est_range[0],
                xmax = est_range[1], start_params = start_params[0:nparams - 2])[0]

        print_res("Results from " + numToWords(nstates - 1) + 
                " state fit ", start_params[0:nparams - 2], level = "INFO")

        # This is necessary to ensure that fit_ansatz_array gives the correct data.
        self.gen_fit_data(xmin, min(int(self._Nt/2), xmax), correlated = False)

        higher_state_data = self._fit_ydata - self.fit_ansatz_array(
                start_params[0:nparams - 2])

        sign_changed = False
        if higher_state_data[0] < 0:
            lg.info("Higher state seems to have negative amplitude!")
            higher_state_data *= -1
            sign_changed = True

        lg.info("Estimate parameters: Linear fit for higher state...")
        tmp_start_params, tmp_start_params_cov = curve_fit(linear, self._fit_xdata,
                np.log(np.abs(higher_state_data)),
                sigma=1 / higher_state_data * self._fit_edata)


        start_params[nparams - 2] = np.exp(tmp_start_params[0]) * \
            np.exp(-tmp_start_params[1] * self._Nt / 2) * 2 / self._mult_const

   
        # We already multiplied with -1, therefore > instead of <
        if sign_changed:
            start_params[nparams - 2] *= -1

        start_params[nparams - 1] = tmp_start_params[1]
     
        lg.info("Estimate parameters: Start", numToWords(nstates),
            " state fit with fixed parameters...")

        def red_corr(x, params, fixed_params, nstates):
            return self._corr(x, fixed_params, nstates - 1, 0)\
                    + self._corr(x, params, 1, 0)


        self._func = red_corr
        self._args = (start_params[0:nparams - 2], nstates)
        self._grad_args = (1,0)
        self._hess_args = (1,0)

        print_res("Starting parameters for first " + numToWords(nstates) + " state fit", 
                start_params, level = "INFO")

        start_params[nparams - 2:], tmp_res_err, chi_dof = self.try_fit(self._std_algs,
                xmin = xmin, xmax = xmax, start_params = start_params[nparams - 2:],
                correlated = False)


        # Reset function
        self._func = self._corr
        self._args = (nstates, 0)
        self._grad_args = (nstates, 0)
        self._hess_args = (nstates, 0)

        return start_params





    """Estimate start parameters for one state fit with oscillations"""
    def _est_params_one_state_osc(self, xmin, xmax):
        even = self._xdata % 2 == 0
        odd = self._xdata % 2 == 1


        try:
            self.gen_fit_data(xmin, xmax, correlated = False, ind = even)

            # Prevent parent class from generating new fit data, as we set them manually here
            self._no_new_data = True

            lg.info("Fit even part...")
            start_params = self._est_params_one_state(xmin, xmax, ind = even)
            res_e = self.simple_corr_fit(xmin, xmax, start_params = start_params,
                    correlated = False, nstates = 1, nstates_osc = 0)[0]
            print_res("Results from even fit", res_e, level = "INFO")


            self.gen_fit_data(xmin, xmax, correlated = False, ind = odd)

            lg.info("Fit odd part...")
            start_params = self._est_params_one_state(xmin, xmax, ind = odd)
            res_o = self.simple_corr_fit(xmin, xmax, start_params = start_params,
                    correlated = False, nstates = 1, nstates_osc = 0)[0]
            print_res("Results from odd fit", res_o, level = "INFO")

            self._no_new_data = False

            # In principle res_e should be larger than res_o, but there are cases where this is
            # not the case. We have make sure that the amplitude is positive at his point
            A_osc = abs(res_e[0] - res_o[0]) / 2
            A_no = (res_o[0] + res_e[0]) / 2

            m = (res_e[1] + res_o[1]) / 2

            start_params = np.array([A_no, m, A_osc, m])
        finally:
            self._no_new_data = False
        return start_params




    """Estimate start parameters for multiple state fit with oscillations"""
    def _est_params_osc(self, xmin, xmax, nstates, nstates_osc):

        lg.info("Estimate parameters: Fit range for %d + %d state fit: ["
                % (nstates, nstates_osc), xmin, ",", xmax, "]")

        est_range = self._get_est_range(xmin, xmax,
                self.est_weights_osc[max(nstates, nstates_osc)])
        nparams_no = 2 * nstates
        nparams = 2 * (nstates + nstates_osc)


        start_params = np.ones(nparams)

        if nstates + nstates_osc == 2:
            return self._est_params_one_state_osc(xmin, xmax)
        else:
            if nstates <= nstates_osc:
                tmp_params = self._est_params_osc(est_range[0], est_range[1], nstates,
                    nstates_osc - 1)

                print_res("Start parameters for %d + %d fit"
                        % (nstates, nstates_osc - 1), tmp_params, level = "INFO")

                tmp_params = self.simple_corr_fit(*est_range, correlated = False,
                        start_params = tmp_params,
                        nstates = nstates, nstates_osc = nstates_osc - 1)[0]

                start_params[:nparams - 2] = tmp_params

                m = 5/4. * start_params[-3]
                A = self._ydata[int(xmin)] - self._func(xmin,
                        tmp_params, nstates, nstates_osc - 1)
                A /= - np.cos(np.pi*xmin) * np.cosh(m * (xmin - self._Nt/2)) * self._mult_const

                start_params[-2] = A
                start_params[-1] = m
            else:
                tmp_params = self._est_params_osc(est_range[0], est_range[1], nstates - 1,
                    nstates_osc)

                print_res("Start parameters for %d + %d fit"
                        % (nstates - 1, nstates_osc), tmp_params, level = "INFO")
                tmp_params = self.simple_corr_fit(*est_range, correlated = False,
                        start_params = tmp_params,
                        nstates = nstates - 1, nstates_osc = nstates_osc)[0]
                

                start_params[:nparams_no - 2] = tmp_params[:nparams_no - 2]
                start_params[nparams_no:] = tmp_params[nparams_no - 2:]

                m = 5/4. * start_params[nparams_no - 3]
                A = self._ydata[int(xmin)] - self._func(xmin,
                        tmp_params, nstates - 1, nstates_osc)
                A /= np.cosh(m * (xmin - self._Nt/2)) * self._mult_const

                start_params[nparams_no - 2] = A
                start_params[nparams_no - 1] = m


        return start_params


    def plot_corr(self, filename, params, params_err, ranges,
                 title=None, notex=False, plot_ylog=True, xmin = 1, norm_func = None,
                 size = (15,10), no_error = True, font_size = 16, **kwargs):

        mult_const = self._mult_const
        self._mult_const = 1


        if np.max(ranges) <= self._Nt/2:
            xmax = self._Nt/2
        else:
            xmax = np.max(self._xdata)


        args_data = {'xmin' : xmin, 'xmax' : xmax, 'title' : title}
        args_data = merge_two_dicts(args_data, kwargs)
        self.plot_fit(filename, params, params_err, notex = notex, ylog = plot_ylog,
                norm_func = norm_func, args_data = args_data, font_size = font_size,
                ranges = ranges, size = size, no_error = no_error, fix_ylim = True)

        self._mult_const = mult_const



    def get_func(self, x, params = None, params_err = None):
        if params is not None:
            params = self._apply_mult_const(params)
        if params_err is not None:
            params_err = self._apply_mult_const(params_err)
        return super().get_func(x, params, params_err)






   
    def _corr(self, x, params, nstates, nstates_osc):
        if len(params) != 2 * (nstates + nstates_osc):
            raise IndexError("_corr: Parameters do not agree with number of states")

        dx = x - self._Nt / 2

        ret = np.zeros_like(x, dtype = float)
        params = np.asarray(params)
        for i in range(0, nstates):
            i *= 2 # because we have two parameters per state
            ret += params[i] * np.cosh(params[i+1] * dx)

        for i in range(nstates, nstates + nstates_osc):
            i *= 2 # because we have two parameters per state
            ret -= np.cos(np.pi*x) * params[i] * np.cosh(params[i+1] * dx)

        return self._mult_const * ret


    def _grad_corr(self, x, params, nstates, nstates_osc):
        if len(params) != 2 * (nstates + nstates_osc):
            raise IndexError("_grad_corr: Parameters do not agree with number of states")

        dx = x - self._Nt / 2

        ret = np.zeros((len(params), len(x)))
        for i in range(0, nstates):
            i *= 2 # we have two parameters per state
            ret[i,:] = np.cosh(params[i+1] * dx)
            ret[i+1,:] = params[i] * np.sinh(params[i+1] * dx) * dx

        for i in range(nstates, nstates + nstates_osc):
            i *= 2 # we have two parameters per state
            ret[i,:] = -np.cos(np.pi*x) * np.cosh(params[i+1] * dx)
            ret[i+1,:] = -np.cos(np.pi*x) * params[i] * np.sinh(params[i+1] * dx) * dx
        return self._mult_const * ret





    def _hess_corr(self, x, params, nstates, nstates_osc):
        if len(params) != 2 * (nstates + nstates_osc):
            raise IndexError("_hess_corr: Parameters do not agree with number of states")

        dx = x - self._Nt / 2

        ret = np.zeros((len(params), len(params), len(x)))

        for i in range(0, nstates):
            i *= 2 # we have two parameters per state
            ret[i,i+1,:] = self._mult_const * np.sinh(params[i+1] * dx) * dx
            ret[i+1,i,:] = ret[i,i+1]
            ret[i+1,i+1,:] = self._mult_const * params[i] * np.cosh(params[i+1] * dx) * dx**2

        for i in range(nstates, nstates + nstates_osc):
            i *= 2 # we have two parameters per state
            ret[i,i+1,:] = -self._mult_const * np.cos(np.pi*x) * np.sinh(params[i+1] * dx) * dx
            ret[i+1,i,:] = ret[i,i+1]
            ret[i+1,i+1,:] = -self._mult_const * np.cos(np.pi*x) * params[i] \
                            * np.cosh(params[i+1] * dx) * dx**2
        return ret



def linear(x, a, b):
    return a - b * x