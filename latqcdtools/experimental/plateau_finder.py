import numpy as np
from latqcdtools.weighted_average import bootstr_add_dist
import matplotlib.pyplot as plt
from latqcdtools.plotting import *
from latqcdtools.autoscale import auto_range_err
from latqcdtools.tools import remove_nan, remove_large, remove_large_err

def find_nearest_upper(array, value):
    if value >= np.max(array):
        return np.max(array)
    array = np.asarray(array)
    array = array[array > value]
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_lower(array, value):
    if value <= np.min(array):
        return np.min(array)
    array = np.asarray(array)
    array = array[array <= value]
    idx = (np.abs(array - value)).argmin()
    return array[idx]


class PlateauFinder:

    def __init__(self, xdata, ydata, edata, amp = None, amp_err = None, chi_dof = None, **kwargs):

        self._acc_factor = kwargs.get('acc_factor', 0.5)
        self._filename = kwargs.get('filename', None)
        self._xdata_col = kwargs.get('xdata_col', 1)
        self._data_col = kwargs.get('data_col', 2)
        self._edata_col = kwargs.get('edata_col', None)
        self._xmin = kwargs.get('xmin', -np.inf)
        self._xmax = kwargs.get('xmax', np.inf)
        self._chi_col = kwargs.get('chi_col', None)
        self._out_name = kwargs.get('out_name', None)
        self._title = kwargs.get('title', None)
        self._npoints = kwargs.get('npoints', None)
        self._show_plot = kwargs.get('show_plot', False)
        self._hist_name = kwargs.get('hist_name', None)
        self._xlabel = kwargs.get('xlabel', None)
        self._ylabel = kwargs.get('ylabel', None)
        self._amp_label = kwargs.get('amp_label', None)
        self._no_tex = kwargs.get('no_tex', False)
        self._amp_col = kwargs.get('amp_col', None)
        self._auto_range = kwargs.get('auto_range', False)
        self._err_threshold = kwargs.get('err_threshold', 0.5)


        rm_list = [xdata, ydata, edata]


        if chi_dof is not None:
            rm_list.append(chi_dof)

        if amp is not None:
            rm_list.append(amp)
        if amp_err is not None:
            rm_list.append(amp_err)


        ind = (xdata < self._xmax) & (xdata > self._xmin)
        for i in range(len(rm_list)):
            rm_list[i] = rm_list[i][ind]

        rm_list = remove_nan(*rm_list, test_cols = (1,2))
        rm_list = remove_large(*rm_list, test_cols = (1,2))


        self._xdata = np.array(rm_list[0])
        self._data = np.array(rm_list[1])
        self._data_err = np.array(rm_list[2])

        if chi_dof is not None:
            self._chi_dof = np.array(rm_list[3])
        else:
            self._chi_dof = None

        if amp is not None:
            self._amp = np.array(rm_list[3 + (self._chi_col is not None)])
            self._amp_err = np.array(rm_list[4 + (self._chi_col is not None)])
        else:
            self._amp = None
            self._amp_col = None


        if self._npoints is None:
            self._npoints = max(math.ceil(len(self._data) / 2), 2)
        

        if len(self._data) == 0:
            raise ValueError("No data points left")


        self._numb_click = 0
        self._ylim = None





    # Not used at the moment
    def _get_min_err(self):
        err = np.inf
        ind = 0
        for i in range(0, len(self._data_err) - self._npoints):
            e = 0
            for j in range(self._npoints):
                e += self._data_err[i + j]
            if e < err:
                err = e
                ind = i
        return [ind, self._npoints + ind]




    def _get_plateau(self):
        self._numb_clicks = 0

        bounds = [np.nan, np.nan]

        def onclick(event):

            self._numb_click += 1
            if self._numb_click == 1:
                bounds[0] = event.xdata
            if self._numb_click == 2:
                bounds[1] = event.xdata

                bounds.sort()
                bounds[0] = find_nearest_upper(self._xdata, bounds[0])
                bounds[1] = find_nearest_lower(self._xdata, bounds[1])
                plt.close()

        fig, ax = plt.subplots()
        cid = fig.canvas.mpl_connect('button_press_event', onclick)

        line = plot_lines(np.array(self._xdata), self._data, self._data_err, ax=ax,
                linewidth = 2, elinewidth = 2, markersize = 5, capsize = 3)

        limits = auto_range_err(self._data, self._data_err, self._acc_factor)

        plt.ylim(limits)
        set_params(xlabel=self._xlabel, ylabel=self._ylabel)
        f = zoom_factory(ax)
        if self._chi_dof is not None and self._amp is None:
            ax2 = ax.twinx()

            limits_chi = auto_range(self._chi_dof, 1.5)
            plot_lines(np.array(self._xdata), self._chi_dof, color=colors[3],
                    label = "$\chi$/d.o.f")
            plot_lines(np.array(self._xdata),
                    np.ones(len(self._xdata)), marker=None, color=colors[3])
            ax2.set_ylim(limits_chi)
            ax2.set_ylabel("$\\chi/\\mathrm{d.o.f}$")
        elif self._amp is not None:
            ax2 = ax.twinx()
            limits_amp = auto_range_err(self._amp, self._amp_err, 3)
            plot_dots(np.array(self._xdata), self._amp, self._amp_err, color=colors[3],
                    label = self._amp_label)
            ax2.set_ylabel(self._amp_label)
            ax2.set_ylim(limits_amp)

        handles, labels = get_legend_handles()
        if self._ylabel is not None:
            handles.append(line)
            labels.append(self._ylabel)
        plt.legend(handles, labels)

        plt.sca(ax)

        # Show the plot and wait for user input
        plt.show()

        
        limits = ax.get_ylim()

        return bounds, limits




    def _auto_plateau(self):
        diff = np.inf
        ind = 0
        for i in range(0, len(self._data) - self._npoints + 1):
            df = 0
            av, err = bootstr_add_dist(self._data[i:i + self._npoints],
                    self._data_err[i:i + self._npoints])
            df = err

            if df < diff:
                diff = df
                ind = i
        return self._xdata[ind], self._xdata[ind + self._npoints - 1]

            


    def _get_bounds(self, bounds = None):
        if bounds is None:
            if not self._auto_range:
                self._bounds, self._ylim = self._get_plateau()
            else:
                self._bounds = self._auto_plateau()
                self._ylim = auto_range_err(self._data, self._data_err, self._acc_factor)

        else:
            bounds = list(bounds)
            bounds.sort()
            bounds[0] = find_nearest_upper(self._xdata, bounds[0])
            bounds[1] = find_nearest_lower(self._xdata, bounds[1])
            self._bounds = bounds
            self._ylim = auto_range_err(self._data, self._data_err, self._acc_factor)




    def get_average(self, bounds = None):
        self._get_bounds(bounds)
        if self._hist_name is not None:
            if self._no_tex:
                init_notex()
            else:
                latexify()

        ind = (self._xdata>= self._bounds[0]) & (self._xdata<= self._bounds[1])
        mean, err = bootstr_add_dist(
                *remove_large_err(self._data[ind], self._data_err[ind], threshold = self._err_threshold,
                    col_val = 0, col_err = 1),
                plot_hist = (self._hist_name is not None))


        if self._hist_name is not None:
            set_params(xlabel = "$m$", title = self._title)
            plt.savefig(self._hist_name)

        self._mean = mean
        self._err = err

        if self._amp_col is not None:
            amp_mean, amp_err = bootstr_add_dist(
                    self._amp[ind],
                    self._amp_err[ind])

            return mean, err, amp_mean, amp_err, self._bounds[0], self._bounds[1]
        else:
            return mean, err, self._bounds[0], self._bounds[1]



    def plot_plat(self):
        if self._out_name is not None:
            if self._no_tex:
                init_notex()
            else:
                latexify()
        ax=plt.gca()
        plot_lines(np.array(self._xdata), self._data, self._data_err, marker='.', markersize=4,
                title = self._title, label=self._ylabel)

        plt.ylim(self._ylim)
        plot_fill(self._bounds,
                [self._mean, self._mean], [self._err, self._err],
                xlabel=self._xlabel, ylabel=self._ylabel, color=colors[0], alpha=0.3)

        if self._chi_dof is not None and self._amp is None:
            ldg = ax.legend(*get_legend_handles(), loc=2, bbox_to_anchor=(0, -.1), borderaxespad=0.)
            ldg.get_frame().set_linewidth(0.5)
            ldg.get_frame().set_alpha(0)
            ldg.get_frame().set_linewidth(0.0)

            clear_legend_labels()

            ax2 = ax.twinx()
            line = plot_lines(np.array(self._xdata), self._chi_dof, color=colors[3], ax=ax2,
                    label="$\\chi/\\mathrm{d.o.f}$")
            xlims = ax.get_xlim()
            plot_lines([xlims[0], xlims[1]], [1, 1], marker=None, color=colors[3],
                    ax=ax2, alpha_lines=0.3)

            limits = auto_range(self._chi_dof, 1.5)
            ax2.set_ylabel("$\\chi/\\mathrm{d.o.f}$")
            ax2.set_ylim(limits)
            ldg = ax2.legend(*get_legend_handles(), bbox_to_anchor=(1.1, -.1), borderaxespad=0., ncol = 2)
            ldg.get_frame().set_alpha(0)
            ldg.get_frame().set_linewidth(0.0)
            ldg.get_frame().set_linewidth(0)
        elif self._amp is not None:
            ldg = ax.legend(*get_legend_handles(), loc=2, bbox_to_anchor=(0, -.1), borderaxespad=0.)
            ldg.get_frame().set_linewidth(0.5)
            ldg.get_frame().set_alpha(0)
            ldg.get_frame().set_linewidth(0.0)

            clear_legend_labels()

            ax2 = ax.twinx()
            line = plot_dots(np.array(self._xdata), self._amp, self._amp_err, color=colors[3], ax=ax2,
                    label=self._amp_label)
            xlims = ax.get_xlim()

            limits = auto_range(self._amp, 3)
            ax2.set_ylabel(self._amp_label)
            ax2.set_ylim(limits)
            ldg = ax2.legend(*get_legend_handles(), bbox_to_anchor=(1.1, -.1), borderaxespad=0., ncol = 2)
            ldg.get_frame().set_alpha(0)
            ldg.get_frame().set_linewidth(0.0)
            ldg.get_frame().set_linewidth(0)
        set_params(title = self._title)

        if self._out_name is not None:
            plt.savefig(self._out_name)
        if self._show_plot:
            plt.show()