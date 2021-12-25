import matplotlib as mpl
import socket
from matplotlib.ticker import AutoMinorLocator
from latqcdtools.tools import check_numpy
from latqcdtools.autoscale import auto_range
from latqcdtools.statistics import error_prop_func, norm_cov
import latqcdtools.logger as logger
import math as math
if (socket.gethostname().startswith("gt") or socket.gethostname().startswith("gx")
        or socket.gethostname().startswith("p0")):
    mpl.use('Agg')


import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
import itertools
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.colors as cl


"""Collection of functions for plotting different data sets using the matplotlib"""


zod = 1 # use for zorder in plot commands will be increased
colors_1 = ['#d32d11', '#0081bf', '#e5af11', '#7c966d', '#7570b3', '#ff934f', '#666666', '#D186B3']
colors_2 = ['#396AB1', '#DA7C30', '#3E9651', '#CC2529', '#535154', '#6B4C9A', '#922428', '#948B3D']
colors_3 = ['#A6CEE3','#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00']

colors = colors_1

markers_1 = ['o', 'v', 'D', 's', 'p', '^', 'h', 'd', 'x', '+', '*']

markers = itertools.cycle(markers_1)

legend_handles = []
legend_labels = []

def set_markers(marker_set = markers_1):
    global markers
    markers = itertools.cycle(marker_set)


def set_colors(color_set = colors_1):
    global colors
    colors = color_set


default_params = {
# Column for the xdata when plotting a file
        'xcol': None,

# Column for the ydata when plotting a file
        'ycol': None,

# Column for the errors in y-direction when plotting a file
        'yecol': None,

# Column for the errors in x-direction when plotting a file
        'xecol': None,

# When passing extra parameter to a function this key defines whether the parameters
# get expanded or not. i.e func(x *param, or func(xparam,
        'expand': True,

# Axis object that is used for the plots. Default is matplotlib.pyplot
        'ax': plt,

# Style when plotting a file
        'style': "dots",

# X-label
        'xlabel': None,

# Y-label
        'ylabel': None,

# Label in legend
        'label': None,

# Title of the plot
        'title': None,

# Transperancy for different plots
        'alpha': 0.5,

# Transperancy for different dots
        'alpha_dots': None,

# Transperancy for different lines
        'alpha_lines': 1,

# Transperancy for edges of error bands
        'alpha_fill_edge': 1,

# Transperancy for labels
        'alpha_label': 0,

# Transperancy for the legend
        'alpha_legend': 0,

# Transperancy for labels
        'alpha_label': 0,

# Linewidth of the legend border
        'linewidth_legend': 0,

# Color for different plots
        'color': None,

# Scale data in xdata by this factor
        'xscale': 1.0,

# Scale data in ydata by this factor
        'yscale': 1.0,

# Number of points for function plotting
        'npoints' : 1000,

# Position of the legend
#  ======================================
#  loc positions:
# 'upper right':  : 1
# 'upper left':   : 2
# 'lower left':   : 3
# 'lower right':  : 4
# 'right':        : 5
# 'center left':  : 6
# 'center right': : 7
# 'lower center': : 8
# 'upper center': : 9
# 'center':       : 10
#  ======================================
        'legendpos': 'best',

# Manual position of the legend
        'bbox_to_anchor': None,

# Number of columns in the legend
        'legend_ncol': 1,

# Spacing between columns in the legend
        'legend_col_spacing': None,

# Spacing between symbol and text in legend
        'handletextpad' : 0.2,

# Linewidth of samples in the legend
        'legend_linewidth': None,

# Title of the legend
        'legend_title': None,
# Shift legend label text
        'shift_legend_label_text_x':None,
        'shift_legend_label_text_y':None,

# Minimum xvalue to be plotted, i.e. xmin is the lower bound on the x-data
# that will be used for plotting. This does NOT directly change the x-range.
        'xmin': None,

# Similarly, maximium xvalue to be plotted
        'xmax': None,

# Marker for plotting dots
        'marker': "iter",

# Size of the dots
        'markersize': 3.5,

# Linewidth of line plots
        'linewidth': 0.1,

# Position of the sub_plot_location
        'loc': 1,

# First edge of the connection line from the zoom window to sub plot
        'loc1': 4,

# Second edge of the connection line from the zoom window to sub plot
        'loc2': 2,

# Padding between subplot and major plot axis
        'borderpad': 0.5,

# Length of caps af error bars
        'capsize': 1.5,

# Linewidth of the error bars of caps af error bars
        'elinewidth': 0.5,

# Put xlabel and ylabel into the plotting area
        'labelsintoplot': False,

# If labelsintplot is true use this to shift the positions of the xlabels
        'xlabelpos': None,

# If labelsintplot is true use this to shift the positions of the xlabel
        'ylabelpos': None,

# Default font size for all kind of text
        'font_size' : 9,
        
# Fill color of points. Set to None (not as string) to have filled symbols
        'point_fill_color': "None",

# z order of plot
        'zod' : None,
        }



def set_default_param(**kwargs):
    for key, val in kwargs.items():
        default_params[key] = val



def fill_param_dict(params):

    """Collection of default parameters for plotting routines. If a key does not exist in 
    params, it is defined with a default value
    
    Parameters
    ----------

    params: dictionary
        Dictionary with all parameters that are already set by the user

    Returns
    -------
        Dictionary filled with all default parameters
    """
    if 'show_leg' not in params:
        # When filling params for the first time, we check if we will show the legend.
        # This is triggered by one of the following keys
        for key in ('legend_title', 'legendpos', 'legend_ncol',
                'legendpos_col_spacing', 'legend_linewidth', 'label', 'alpha_legend'):
            if key in params:
                params['show_leg'] = True
    if 'show_leg' not in params:
        params['show_leg'] = False

    for key, val in default_params.items():
        params.setdefault(key,val)


    



def add_optional(params):
    """Optional parameter that are not defined in fill_param_dict are collected by this
    function
    Parameters
    ----------

    params: dictionary
    Dictionary with parameters set by the user

    Returns
    -------
    Dictionary with optional parameters
    """
    reference = {}
    fill_param_dict(reference)
    reference['show_leg'] = False
    ret = {}
    for key in params:
        if not key in reference:
            ret[key] = params[key]
    return ret


def set_mxtics(n):
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(n))

def set_mytics(n):
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(n))

# Try xycoords = 'data' for data coordinates
def set_label(label, pos, color="black", xycoords = 'axes fraction', zorder = None, **params):

    fill_param_dict(params)
    optional = add_optional(params)

    font_size = params['font_size']

    plt.gca().annotate(label, xy=pos, xycoords=xycoords, color = color,
            fontsize = font_size, bbox=dict(linewidth = 0, facecolor = 'white',
                alpha=params['alpha_label']), zorder = zorder, **optional)



def set_params(ax = plt, **params):

    """
    Set additional parameters to the plot. For example set a title or label
    Parameters
    ----------
    ax : Axis object, optional, default = plt
    The axis to which the parameters should be set
    **params :
    Additional parameters that can be set.
    
    """

    fill_param_dict(params)

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    if params['xlabel'] is not None:
        if params['labelsintoplot'] or params['xlabelpos'] is not None:
            if params['xlabelpos'] is None:
                params['xlabelpos'] = (0.95,0.027)
            set_label(params['xlabel'], params['xlabelpos'], ha='right',
                    va='bottom', xycoords='axes fraction', alpha_label = params['alpha_label'],
                    font_size = params['font_size'], zorder = zod)
        else:
            if ax == plt:
                ax.xlabel(params['xlabel'])
            else:
                ax.set_xlabel(params['xlabel'])

    if params['ylabel'] is not None:
        if params['labelsintoplot'] or params['ylabelpos'] is not None:
            if params['ylabelpos'] is None:
                params['ylabelpos'] = (0.025,0.972)
            set_label(params['ylabel'], params['ylabelpos'], ha='left', va='top',
                    xycoords='axes fraction', alpha_label = params['alpha_label'],
                    font_size = params['font_size'], zorder = zod)
        else:
            if ax == plt:
                ax.ylabel(params['ylabel'])
            else:
                ax.set_ylabel(params['ylabel'])

    if params['title'] is not None:
        if ax == plt:
            ax.title(params['title'])

    if params['show_leg']:

        leg = ax.legend(legend_handles, legend_labels, numpoints=1, bbox_to_anchor = params['bbox_to_anchor'],
                title=params['legend_title'], loc=params['legendpos'],
                ncol=params['legend_ncol'], columnspacing=params['legend_col_spacing'],
                handletextpad = params['handletextpad'])

        leg.get_frame().set_alpha(params['alpha_legend'])
        leg.get_frame().set_linewidth(params['linewidth_legend'])
        leg.set_zorder(zod)
        if params['legend_linewidth'] is not None:
            for legobj in leg.legendHandles:
                legobj.set_linewidth(params['legend_linewidth'])

        if params['shift_legend_label_text_y']: # or params['shift_legend_label_text_y'] :
            x_shift_ = params['shift_legend_label_text_x'] 
            y_shift_ = params['shift_legend_label_text_y'] 
            if x_shift_ is None:
                x_shift_ = 0
            if y_shift_ is None:
                y_shift_ = 0
            renderer = plt.gcf().canvas.get_renderer()

            for i in range(len(legend_handles)):
              x_shift = x_shift_*leg.texts[i].get_window_extent(renderer).width
              y_shift = y_shift_*leg.texts[i].get_window_extent(renderer).height
              leg.texts[i].set_position((x_shift, y_shift))


def plot_file(filename, xcol=1, ycol=2, yecol=None, xecol=None, ax=plt, func = None,
        func_args = (), **params):
    fill_param_dict(params)
    data = np.loadtxt(filename, dtype = np.str).transpose()
    if xcol is not None:
        xdata = np.array(data[xcol - 1], dtype = float)
    else:
        xdata = np.arange(len(ydata))
    ydata = np.array(data[ycol - 1], dtype = float)
    yedata = None
    xedata = None

    if yecol is not None:
        yedata = np.array(data[yecol - 1], dtype = float)
    else:
        if params['style'] == "fill":
            raise ValueError("Need error column for filled plotting")
    if xecol is not None:
        xedata = np.array(data[xecol - 1], dtype = float)


    if func is not None:
        if yedata is not None:
            if xedata is not None:
                xdata, ydata, yedata, xedata = func(xdata, ydata, yedata, xedata, *func_args)
            else:
                xdata, ydata, yedata = func(xdata, ydata, yedata, *func_args)
        else:
            xdata, ydata = func(xdata, ydata, *func_args)

    if params['style'] == "dots":
        return plot_dots(xdata, ydata, yedata=yedata, xedata=xedata, **params)
    if params['style'] == "lines":
        return plot_lines(xdata, ydata, yedata=yedata, xedata=xedata, **params)
    if params['style'] == "fill":
        return plot_fill(xdata, ydata, yedata=yedata, **params)
    if params['style'] == "band":
        return plot_band(xdata, ydata, yedata, xedata, **params)

    raise ValueError("No such style: " + params['style'])


def set_auto_range(scale_factor = 0.3):
    data_list = []

    for line in plt.gca().get_lines():
        data_list += list(line.get_ydata())

    ymin, ymax = auto_range(data_list, scale_factor)
    set_ymin(ymin)
    set_ymax(ymax)



def plot_dots(xdata, ydata, yedata = None, xedata = None, **params):
    xdata, ydata, yedata, xedata = check_numpy(xdata, ydata, yedata, xedata)
    fill_param_dict(params)
    optional = add_optional(params)

    xdata, ydata, yedata, xedata = remove_points(xdata, ydata, yedata, xedata,
            xmin = params['xmin'], xmax = params['xmax'])

    ax = params['ax']
    xscale = params['xscale']
    yscale = params['yscale']
    if xedata is not None:
        xedata=np.copy(xedata*xscale)
    if yedata is not None:
        yedata=np.copy(yedata*yscale)

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    marker = params['marker']
    if marker == "iter":
        marker = next(markers)
    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale,
                yerr=yedata, xerr=xedata, marker=marker,
                linestyle='None', linewidth=params['linewidth'],
                alpha=params['alpha_dots'],  color=params['color'], zorder=zod,
                markersize=params['markersize'], capsize=params['capsize'],
                elinewidth=params['elinewidth'], markeredgewidth=params['elinewidth'],
                markerfacecolor = params['point_fill_color'], **optional)
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata,
                xerr=xedata, marker=marker, linestyle='None', 
                linewidth=params['linewidth'], alpha=params['alpha_dots'],
                zorder=zod, markersize=params['markersize'], capsize=params['capsize'],
                elinewidth=params['elinewidth'], markeredgewidth=params['elinewidth'],
                markerfacecolor = params['point_fill_color'], **optional)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append(ebar)

    globals()['zod'] += 1
    set_params(**params)
    
    return ebar

def plot_bar(xdata, ydata, width=None, align='edge', alpha=1, edgecolor='#666677',linewidth=0.2, **params):
    if width is None:
        width = xdata[1] - xdata[0]

    if alpha is None:
        alpha=params['alpha']

    xdata, ydata = check_numpy(xdata, ydata)
    fill_param_dict(params)
    optional = add_optional(params)
    ax = params['ax']
    params['alpha'] = 1
    params['linewidth'] = 0

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    if params['color'] is not None:
        bar = ax.bar(xdata, ydata, color=params['color'], zorder=zod, width=width, align=align, edgecolor=edgecolor,
                linewidth=linewidth, alpha=alpha, **optional)
    else:
        bar = ax.bar(xdata, ydata, zorder=zod, width=width, align=align, edgecolor=edgecolor,
                linewidth=linewidth, alpha=alpha, **optional)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append(bar)
        params['shift_legend_label_text_y'] = -0.22
        params['shift_legend_label_text_x'] = 0.05


    globals()['zod'] += 1
    set_params(**params)

    return bar


def plot_hist(data, logx = False, bins = None):
    if logx:
        xmin = np.min(data)
        xmax = np.max(data)

        if bins is None:
            bins = 50
        else:
            bins = bins
        bins = np.logspace(np.log10(xmin),np.log10(xmax), bins)
        plt.hist(data, bins = bins)
        plt.xscale('log')
    else:
        if bins is None:
            bins = 'auto'
        plt.hist(data, bins = bins)



def plot_lines(xdata, ydata, yedata=None, xedata=None, **params):
    xdata, ydata, yedata, xedata = check_numpy(xdata, ydata, yedata, xedata)
    fill_param_dict(params)
    optional = add_optional(params)
    xdata, ydata, yedata, xedata = remove_points(xdata, ydata, yedata, xedata,
            xmin = params['xmin'], xmax = params['xmax'])
    xscale=params['xscale']
    yscale=params['yscale']
    ax = params['ax']

    if xedata is not None:
        xedata=np.copy(xedata*xscale)
    if yedata is not None:
        yedata=np.copy(yedata*yscale)
    ebar = None

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    marker = params['marker']
    if marker == "iter":
        marker = next(markers)



    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale,
                yerr=yedata, xerr=xedata, marker=marker,
                linestyle='None', linewidth=params['linewidth'], color=params['color'],
                zorder=zod, markersize=params['markersize'], capsize=params['capsize'],
                elinewidth=params['elinewidth'], markeredgewidth=params['elinewidth'],
                alpha=params['alpha_dots'], markerfacecolor = params['point_fill_color'])
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale,
                yerr=yedata, xerr=xedata, marker=marker,
                linestyle='None', linewidth=params['linewidth'], 
                zorder=zod, markersize=params['markersize'], capsize=params['capsize'],
                elinewidth=params['elinewidth'], markeredgewidth=params['elinewidth'],
                alpha=params['alpha_dots'], markerfacecolor = params['point_fill_color'])

    globals()['zod'] += 1

    col = ebar[0].get_color()

    line = ax.errorbar(xdata*xscale, ydata*yscale, color = col,
            linewidth=params['linewidth'], zorder = zod, alpha = params["alpha_lines"], **optional)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append((line, ebar))

    set_params(**params)
    return ebar





def save_func(func, filename, args=(), func_err=None, args_err=(), grad = None, func_sup_numpy = False,
        **params):

    fill_param_dict(params)

    xmin = params['xmin']
    xmax = params['xmax']


    if params['expand']:
        wrap_func = lambda x, *args: func(x, *args)
        wrap_func_err = lambda x, *args_err: func_err(x, *args_err)
        wrap_grad = lambda x, *args: grad(x, *args)
    else:
        wrap_func = lambda x, *args: func(x, args)
        wrap_func_err = lambda x, *args_err: func_err(x, args_err)
        wrap_grad = lambda x, *args: grad(x, args)


    if xmin is None:
        for line in plt.gca().lines:
            xmin_new = np.min(line.get_xdata())
            if xmin is None:
                xmin = xmin_new
            if xmin_new < xmin:
                xmin = xmin_new
    if xmax is None:
        for line in plt.gca().lines:
            xmax_new = np.max(line.get_xdata())
            if xmax is None:
                xmax = xmax_new
            if xmax_new > xmax:
                xmax = xmax_new

    if xmin is None:
        xmin = -10
    if xmax is None:
        xmax = 10

    xdata = np.arange(xmin, xmax, (xmax - xmin) / params['npoints'])

    if func_sup_numpy:
        ydata = wrap_func(xdata, *args)
    else:
        ydata = np.array([wrap_func(x, *args) for x in xdata])

    if func_err is not None:
        if func_sup_numpy:
            ydata_err = wrap_func_err(xdata, *args_err)
        else:
            ydata_err = np.array([wrap_func_err(x, *args_err) for x in xdata])

        with open(filename, "w") as fout:
            for i in range(len(xdata)):
                print(xdata[i], ydata[i], ydata_err[i], file = fout)

    elif len(args_err) > 0:
        if grad is None:
            logger.warn("Used numerical derivative!")
            wrap_grad = None

        # Arguments that are part of the error propagation
        tmp_args = tuple(args)[0:len(args_err)]

        # Optional arguments, that are constant and, therefore,
        # not part of the error propagation
        tmp_opt = tuple(args)[len(args_err):]

        if func_sup_numpy:
            ydata_err = error_prop_func(xdata, wrap_func, tmp_args, args_err,
                grad = wrap_grad, args = tmp_opt)
        else:
            ydata_err = np.array([error_prop_func(x, wrap_func, tmp_args, args_err,
                grad = wrap_grad, args = tmp_opt) for x in xdata])

        with open(filename, "w") as fout:
            for i in range(len(xdata)):
                print(xdata[i], ydata[i], ydata_err[i], file = fout)

    else:
        with open(filename, "w") as fout:
            for i in range(len(xdata)):
                print(xdata[i], ydata[i], file = fout)










# To plot an error band with an explicit error function, use func_err. args_err are all parameters for func_err
# To use a numerical derivative, just pass the errors of args to args_err
def plot_func(func, args=(), func_err=None, args_err=(), grad = None, func_sup_numpy = False,
        **params):

    fill_param_dict(params)

    params['marker'] = None
    xmin = params['xmin']
    xmax = params['xmax']
    ax = params['ax']


    if params['expand']:
        wrap_func = lambda x, *args: func(x, *args)
        wrap_func_err = lambda x, *args_err: func_err(x, *args_err)
        wrap_grad = lambda x, *args: grad(x, *args)
    else:
        wrap_func = lambda x, *args: func(x, args)
        wrap_func_err = lambda x, *args_err: func_err(x, args_err)
        wrap_grad = lambda x, *args: grad(x, args)


    if xmin is None:
        for line in plt.gca().lines:
            xmin_new = np.min(line.get_xdata())
            if xmin is None:
                xmin = xmin_new
            if xmin_new < xmin:
                xmin = xmin_new
    if xmax is None:
        for line in plt.gca().lines:
            xmax_new = np.max(line.get_xdata())
            if xmax is None:
                xmax = xmax_new
            if xmax_new > xmax:
                xmax = xmax_new

    if xmin is None:
        xmin = -10
    if xmax is None:
        xmax = 10

    xdata = np.arange(xmin, xmax, (xmax - xmin) / params['npoints'])

    if func_sup_numpy:
        ydata = wrap_func(xdata, *args)
    else:
        ydata = np.array([wrap_func(x, *args) for x in xdata])

    if func_err is not None:
        if func_sup_numpy:
            ydata_err = wrap_func_err(xdata, *args_err)
        else:
            ydata_err = np.array([wrap_func_err(x, *args_err) for x in xdata])

        return plot_fill(xdata, ydata, ydata_err, **params)

    elif len(args_err) > 0:
        if grad is None:
            logger.warn("Used numerical derivative!")
            wrap_grad = None

        # Arguments that are part of the error propagation
        tmp_args = tuple(args)[0:len(args_err)]

        # Optional arguments, that are constant and, therefore,
        # not part of the error propagation
        tmp_opt = tuple(args)[len(args_err):]

        if func_sup_numpy:
            ydata_err = error_prop_func(xdata, wrap_func, tmp_args, args_err,
                grad = wrap_grad, args = tmp_opt)
        else:
            ydata_err = np.array([error_prop_func(x, wrap_func, tmp_args, args_err,
                grad = wrap_grad, args = tmp_opt) for x in xdata])

        return plot_fill(xdata, ydata, ydata_err, **params)
    else:
        return plot_lines(xdata, ydata, yedata=None, xedata=None, **params)








def remove_points(xdata, *args, xmin = -np.inf, xmax = np.inf):
    if xmin is None:
        xmin = -np.inf

    if xmax is None:
        xmax = np.inf

    ind = (xdata>=xmin) & (xdata<=xmax)
    ret = [xdata[ind]]
    for i in args:
        try:
            if i is not None:
                ret.append(i[ind])
            else:
                ret.append(None)
        except (IndexError, TypeError):
            ret.append(i)
    return ret


def plot_fill(xdata, ydata, yedata, **params):
    xdata, ydata, yedata = check_numpy(xdata, ydata, yedata)
    fill_param_dict(params)
    optional = add_optional(params)
    xdata, ydata, yedata = remove_points(xdata, ydata, yedata,
            xmin = params['xmin'], xmax = params['xmax'])
    xscale=params['xscale']
    yscale=params['yscale']
    ax = params['ax']

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, 
                linewidth=params['linewidth'], color=params['color'], zorder=zod,
                alpha = params['alpha_lines'], **optional)
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, 
                linewidth=params['linewidth'], zorder=zod,
                alpha = params['alpha_lines'], **optional)

    globals()['zod'] += 1

    col = ebar[0].get_color()
    pl = ax.fill_between(xdata*xscale, (np.asarray(ydata*yscale) - np.asarray(yedata*yscale)),
            (np.asarray(ydata*yscale) + np.asarray(yedata*yscale)),
            facecolor=col, alpha=params['alpha'], linewidth=0, zorder=1)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append((ebar, pl))

    set_params(**params)
    return ebar,pl

def plot_band(xdata, low_lim, up_lim, center = None, **params):
    fill_param_dict(params)
    optional = add_optional(params)
    if center is not None:
        xdata, low_lim, up_lim, center = remove_points(xdata, low_lim, up_lim, center,
                xmin = params['xmin'], xmax = params['xmax'])
    else:
        xdata, low_lim, up_lim = remove_points(xdata, low_lim, up_lim, 
                xmin = params['xmin'], xmax = params['xmax'])

    xscale=params['xscale']
    yscale=params['yscale']
    ax = params['ax']

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    globals()['zod'] += 1

    if params['color'] is None:
        pl = ax.fill_between(xdata*xscale, yscale*low_lim,
                yscale*up_lim,
                alpha=params['alpha'], linewidth=0, zorder=1)
    else:
        pl = ax.fill_between(xdata*xscale, yscale*low_lim,
                yscale*up_lim, facecolor=params['color'],
                alpha=params['alpha'], linewidth=0, zorder=1)

    col = cl.rgb2hex(pl.get_facecolor()[0])

    
    if params['alpha_lines'] != 0:
        ax.errorbar(xdata*xscale, yscale*low_lim, color = col,
            linewidth=params['linewidth'], zorder = zod, alpha = params["alpha_fill_edge"])
        ax.errorbar(xdata*xscale, yscale*up_lim, color = col,
            linewidth=params['linewidth'], zorder = zod, alpha = params["alpha_fill_edge"])

    ebar = None
    if center is not None:
        ebar = ax.errorbar(xdata*xscale, center*yscale, 
                linewidth=params['linewidth'], color=col, zorder=zod,
                alpha = params['alpha_lines'], **optional)


    if params['label'] is not None:
        legend_labels.append(params['label'])
        if ebar is not None:
            legend_handles.append((ebar, pl))
        else:
            legend_handles.append(pl)

    set_params(**params)
    if ebar is not None:
        return ebar,pl
    else:
        return pl

def set_xlimit(x_limit=None):
    if x_limit is not None:
        ax = plt.gca()
        ax.set_xlim(x_limit)

def set_ylimit(y_limit=None):
    if y_limit is not None:
        ax = plt.gca()
        ax.set_ylim(y_limit)

def set_xmin(x_min=None):
    if x_min is not None:
        ax = plt.gca()
        x1, x2 = ax.get_xlim()
        ax.set_xlim([x_min,x2])


def set_xmax(x_max=None):
    if x_max is not None:
        ax = plt.gca()
        x1, x2 = ax.get_xlim()
        ax.set_xlim([x1,x_max])

def set_ymin(y_min=None):
    if y_min is not None:
        ax = plt.gca()
        y1, y2 = ax.get_ylim()
        ax.set_ylim([y_min,y2])

def set_ymax(y_max=None):
    if y_max is not None:
        ax = plt.gca()
        y1, y2 = ax.get_ylim()
        ax.set_ylim([y1,y_max])



def set_xrange(xmin=None, xmax=None):
    set_xmin(xmin)
    set_xmax(xmax)

def set_yrange(ymin=None, ymax=None):
    set_ymin(ymin)
    set_ymax(ymax)



def get_xrange():
    ax = plt.gca()
    x1, x2 = ax.get_xlim()
    return x1, x2

def get_yrange():
    ax = plt.gca()
    y1, y2 = ax.get_ylim()
    return y1, y2



def init_notex(fig_width=10, fig_height=7, font_size=9, static_margins=False, init=True):
    
    clear_legend_labels()
    plt.close("all")

    left_margin  = 2 / fig_width
    right_margin = 2 / fig_width
    bottom_margin = 1.5 / fig_height
    top_margin = 1.5 / fig_height

    fig_width /= 2.54
    fig_height /= 2.54


    if font_size is None:
        font_size = default_params['font_size']
    
    set_markers()
    plt.rcParams['legend.handlelength'] = 1.5
    if init:
        plt.rcParams['figure.figsize'] = [fig_width, fig_height]
        plt.rcParams['figure.autolayout'] = True
        plt.rcParams['axes.titlesize'] = font_size
        plt.rcParams['savefig.bbox'] = 'standard'
        plt.rcParams['ytick.labelsize'] = font_size
        plt.rcParams['font.size'] = font_size
        plt.rcParams['axes.labelsize'] = font_size
        plt.rcParams['legend.fontsize'] = font_size
    plt.rcParams['xtick.labelsize'] = font_size
    plt.rc('axes', linewidth=0.5)
    if init:
        fig, ax = plt.subplots()
    else:
        fig = plt.gcf()
        ax = plt.gca()
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    x_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
   
    if static_margins:
        # dimensions are calculated relative to the figure size
        x = left_margin    # horiz. position of bottom-left corner
        y = bottom_margin  # vert. position of bottom-left corner
        w = 1 - (left_margin + right_margin) # width of axes
        h = 1 - (bottom_margin + top_margin) # height of axes

        ax.set_positon([x, y, w, h])

    ax.yaxis.set_major_formatter(y_formatter)
    ax.xaxis.set_major_formatter(x_formatter)
    plt.ticklabel_format(style='sci', scilimits=(-3, 4))
    ax.set_prop_cycle(cycler('color', colors))
    return fig, ax


def clear_legend_labels():

    global legend_labels
    legend_labels = []

    global legend_handles
    legend_handles = []

def get_legend_handles():
    return legend_handles, legend_labels
    
def set_legend_handles(legend_handles, legend_labels):
    globals()['legend_handles'] = legend_handles
    globals()['legend_labels'] = legend_labels

def append_to_legend_handles(legend_handle, legend_label):
    legend_handles.append(legend_handle)
    legend_labels.append(legend_label)

def latexify(fig_width=10, fig_height=7, font_size=9, static_margins=False, init=True, init_fig = None,
        auto_layout = True):
    clear_legend_labels()
    if init_fig is None:
        init_fig = init
    if init_fig:
        plt.close("all")
    
    left_margin  = 2 / fig_width
    right_margin = 2 / fig_width
    bottom_margin = 1.5 / fig_height
    top_margin = 1.5 / fig_height

    fig_width /= 2.54
    fig_height /= 2.54

    if font_size is None:
        font_size = default_params['font_size']
    
    set_markers()
    plt.rcParams['legend.handlelength'] = 1.5
    plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
    plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{braket}"
    plt.rcParams['text.usetex'] = True
    if init:
        plt.rcParams['figure.figsize'] = [fig_width, fig_height]
        plt.rcParams['figure.autolayout'] = auto_layout
        plt.rcParams['axes.titlesize'] = font_size
        plt.rcParams['savefig.bbox'] = 'standard'
        plt.rcParams['ytick.labelsize'] = font_size
        plt.rcParams['font.size'] = font_size
        plt.rcParams['axes.labelsize'] = font_size
        plt.rcParams['legend.fontsize'] = font_size
    plt.rcParams['xtick.labelsize'] = font_size
    plt.rc('axes', linewidth=0.5)
    if init_fig:
        fig = plt.figure()
        ax  = fig.add_subplot(111)
    else:
        fig = plt.gcf()
        ax = plt.gca()
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    x_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
   
    if static_margins:
        ax.set_axis_off()
        # dimensions are calculated relative to the figure size
        x = left_margin    # horiz. position of bottom-left corner
        y = bottom_margin  # vert. position of bottom-left corner
        w = 1 - (left_margin + right_margin) # width of axes
        h = 1 - (bottom_margin + top_margin) # height of axes

        ax = fig.add_axes([x, y, w, h])

    ax.yaxis.set_major_formatter(y_formatter)
    ax.xaxis.set_major_formatter(x_formatter)
    plt.ticklabel_format(style='sci', scilimits=(-3, 4))
    ax.set_prop_cycle(cycler('color', colors))
    return fig, ax



def plot_cov(cov, filename = None, title=None, notex=False, ignore_first = 0, norm = True,
        xrange = None, yrange = None, xmin = None, xmax = None, ymin = None, ymax = None,
        xlabel = "$n_{\\tau/\\sigma}$", ylabel = "$n_{\\tau/\\sigma}$"):
    if not notex:
        latexify()
    else:
        init_notex()
    if norm:
        ncov = norm_cov(cov)
    else:
        ncov = cov
    if title is not None:
        plt.title(title)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xmin is None:
        off_x = 0
    else:
        off_x = xmin

    if ymin is None:
        off_y = 0
    else:
        off_y = ymin

    if xrange is None:
        xrange = np.arange(off_x + ignore_first, off_x + len(cov)+1)
    if yrange is None:
        yrange = np.arange(off_y + ignore_first, off_y + len(cov)+1)

    plt.pcolormesh(xrange, yrange, ncov[ignore_first:, ignore_first:], cmap = "Blues")

    if xmin is not None:
        set_xmin(xmin)
    if ymin is not None:
        set_ymin(ymin)
    if xmax is not None:
        set_xmax(xmax)
    if ymax is not None:
        set_ymax(ymax)

    plt.gca().invert_yaxis()

    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=12)
    if filename is not None:
        plt.savefig(filename)



def plot_eig(cov, filename, title=None, notex=False):
    if not notex:
        latexify()
    else:
        init_notex()
    v, w = np.linalg.eig(cov)
    eig_real = np.real(np.sort(v))
    eig_imag = np.imag(np.sort(v))
    plt.yscale('log')
    plot_bar(range(len(eig_real), 0, -1), eig_real, color='#d32d11', label="real",
            alpha=0.7, title=title, xlabel = "$i$", ylabel = "$E_i$")
    if np.min(eig_imag) != 0:
        plot_bar(range(len(eig_imag), 0, -1), eig_imag, color='#0081bf', label="imag",
                alpha=0.7, title=title, xlabel = "$i$", ylabel = "$E_i$" )
    plt.savefig(filename)
    plt.clf()


# if autoscale is on, call these functions AFTER plt.plot(...) !


def plot_horiz(y, ax = plt, **params):
    fill_param_dict(params)
    color = params['color']
    optional = add_optional(params)
    if color is None:
        color = "black"
    ax.axhline(y = y, linewidth = params['linewidth'], color = color, **optional)

def plot_vert(x, **params):
    fill_param_dict(params)
    optional = add_optional(params)
    color = params['color']
    if color is None:
        color = "black"
    plt.axvline(x = x, linewidth = params['linewidth'], color = color, **optional)

def rel_x(x):
    x1, x2 = plt.gca().get_xlim()
    return x1 + ((x2 - x1) * x)


def rel_y(y):
    y1, y2 = plt.gca().get_ylim()
    return y1 + ((y2 - y1) * y)

def computeTicks (x, step = 5):
    """
    Computes domain with given step encompassing series x
    Parameters
    x: array_like
        A list-like object of integers or floats
    step: optional
        Tick frequency
    """
    xMax, xMin = math.ceil(max(x)), math.floor(min(x))
    dMax, dMin = (xMax + abs((xMax % step) - step) + (step if (xMax % step != 0) else 0),
            xMin - abs((xMin % step)))
    return range(dMin, dMax, step)




# from http://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel

def zoom_factory(ax, base_scale=2.):
    def zoom_fun(event):
        # get the current x and y limits
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0]) * .5
        cur_yrange = (cur_ylim[1] - cur_ylim[0]) * .5
        inv = ax.transData.inverted()
        xdata, ydata = inv.transform((event.x, event.y))
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1 / base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            logger.warn(event.button)
        # set new limits
        ax.set_ylim([ydata - cur_yrange * scale_factor,
                     ydata + cur_yrange * scale_factor])
        plt.draw()  # force re-draw

    fig = ax.get_figure()  # get the figure of interest
    # attach the call back
    fig.canvas.mpl_connect('scroll_event', zoom_fun)

    # return the function
    return zoom_fun


# if you want to plot a "zoom-plot" inside your plot, get this axis object
def zoom_axis(ax, width, height, zx_min, zx_max, zy_min, zy_max, loc=1, loc1=2, loc2=4, borderpad=0.5):
    # width and height with respect to the parent axis should be passed like: width="70%" etc...
    # loc=1 == upper-right
    axins = inset_axes(ax, width, height, loc=loc, borderpad=borderpad) 
    axins.set_xlim(zx_min, zx_max) # apply the x-limits
    axins.set_ylim(zy_min, zy_max) # apply the y-limits

    global zod
    mark_inset(ax, axins, loc1=loc1, loc2=loc2, fc="none", ec="0.5", linewidth=0.5, zorder=zod)
    zod+=1

    return axins # use this axis to plot inside the box!


# plot zoom with given data
def plot_data_zoom(width, height, zx_min, zx_max, xdata, ydata,
        yedata=None, xedata=None, zy_min=None, zy_max=None, **params):
    fill_param_dict(params)
    xscale=params['xscale']
    yscale=params['yscale']
    if params['ax'] is plt:
        params['ax'] = plt.gca()


    if zy_max is None:
        if yedata is None:
            tmp = np.copy(ydata)
        else:
            tmp = ydata+yedata
        zy_max = np.min(tmp)
        for i in range(len(xdata)):
            if zx_max > xdata[i] > zx_min:
                if zy_max < tmp[i]:
                    zy_max = tmp[i]

    if zy_min is None:
        if yedata is None:
            tmp = np.copy(ydata)
        else:
            tmp = ydata-yedata
        zy_min = np.max(tmp)
        for i in range(len(xdata)):
            if zx_max > xdata[i] > zx_min:
                if zy_min > tmp[i]:
                    zy_min = tmp[i]


    

    params['ax'] = zoom_axis(params['ax'], width, height, zx_min*xscale,
            zx_max*xscale, zy_min*yscale, zy_max*yscale,
            loc = params['loc'], loc1 = params['loc1'],
            loc2 = params['loc2'], borderpad=params['borderpad'])

    if params['style'] == "dots":
        return plot_dots(xdata, ydata, yedata=yedata, xedata=xedata, **params), params['ax']
    if params['style'] == "lines":
        return plot_lines(xdata, ydata, yedata=yedata, xedata=xedata, **params), params['ax']
    if params['style'] == "fill":
        return plot_fill(xdata, ydata, yedata=yedata, **params), params['ax']


# plot zoom with given file
def plot_file_zoom(width, height, zx_min, zx_max, filename, xcol=None,
        ycol=None, yecol=None, xecol=None, zy_min=None, zy_max=None, **params):
    fill_param_dict(params)
    xscale=params['xscale']
    yscale=params['yscale']

    data = np.loadtxt(filename).transpose()



    if xcol is not None:
        xdata = data[xcol - 1]
    else:
        xdata = np.arange(len(ydata))
    if ycol is not None:
        ydata = data[ycol - 1]
    else:
        ydata = np.arange(len(xdata))
    yedata = None
    xedata = None
    if yecol is not None:
        yedata = data[yecol - 1]
    if xecol is not None:
        xedata = data[xecol - 1]

    return plot_data_zoom(width, height, zx_min, zx_max, xdata, ydata,
            yedata, xedata, zy_min, zy_max, **params)




class careful_autoscale:
    def __init__(self, ax=None):
        self.ax = ax if ax else plt.gca()
    
    def __enter__(self):
        self.xl = self.ax.get_xlim()
        self.yl = self.ax.get_ylim()
        self.lines = self.ax.get_lines()
        self.lines_visibility = [ l.get_visible() for l in self.lines ]
        [ l.set_visible(False) for l in self.lines ]
        self.ax.autoscale(True)

    def __exit__(self, type, value, traceback):
        self.ax.relim(visible_only=True)
        self.ax.autoscale_view()
        [ l.set_visible(v) for l,v in zip(self.lines, self.lines_visibility) ]