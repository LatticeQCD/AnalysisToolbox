#
# plotting.py
#
# H. Sandmeyer, D. Clarke
#
# Collection of convenience tools for plotting using matplotlib.
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import itertools
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import error_prop_func, norm_cov

# TODO: put in some docstrings for stuff that doesn't have it yet.
#       finally check all the gcas, they may not be needed anymore. also you need to like make sure that all the
#       external methods here have some sort of test; you changed a lot of things and it's important to verify that
#       the changes you made are stable. also expand a bit in the documentation


zod        = 1
initialize = True


colors_1 = ['#d32d11', '#0081bf', '#e5af11', '#7c966d', '#7570b3', '#ff934f', '#666666', '#D186B3']
colors_2 = ['#396AB1', '#DA7C30', '#3E9651', '#CC2529', '#535154', '#6B4C9A', '#922428', '#948B3D']
colors_3 = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00']
colors   = colors_1


markers_1 = ['o', 'v', 'D', 's', 'p', '^', 'h', 'd', 'x', '+', '*']
markers   = itertools.cycle(markers_1)


legend_handles = []
legend_labels = []


default_params = {

    # Hidden python options. (Generally speaking do not touch these.)
    'ax': plt,      # Axis object that is used for the plots. Default is matplotlib.pyplot.
    'expand': True, # Defines whether the parameters get expanded or not, i.e func(x *param) or func(x param).

    # Basic options affecting most plots.
    'xlabel': None,
    'ylabel': None,
    'title': None,
    'label': None,               # What are the data called? (Will appear in legend.)
    'color': None,               # Color for your data. (By default each new set automatically gets different color.)
    'marker': "iter",            # Symbol used for plotting data. (Set to 'None' if you don't want any.)
    'markersize': 3.5,           # Size of the symbols.
    'font_size': 12,             # Default font size for text.
    'font_weight' : 'normal',    # Default style of font ('normal', 'bold', 'heavy', 'light')
    'alpha': 0.5,                # General transparency for data.
    'xscale': 1.0,               # Scale data in xdata by this factor.
    'yscale': 1.0,               # Scale data in ydata by this factor.
    'ticksintoplot' : True,      # Put ticks into plotting area.
    'surroundWithTicks' : True,  # Put ticks also on top and right.
    'labelsintoplot': True,      # Put xlabel and ylabel into plotting area.
    'xlabelpos': None,           # If labelsintplot=True, shift the position (x,y) of the x-label, expressed as percent.
    'ylabelpos': None,
    'alpha_dots': None,          # Transperancy for different dots
    'alpha_lines': 1,            # Transperancy for different lines
    'alpha_fill_edge': 0,        # Transperancy for edges of error bands
    'alpha_label': 1,            # Transperancy for labels
    'linewidth': 1,              # Linewidth of line plots
    'capsize': 1.5,              # Length of caps af error bars
    'elinewidth': 0.5,           # Linewidth of the error bars of caps af error bars
    'point_fill_color': "None",  # Fill color of points. Set to None (not as string) to have filled symbols
    'zod': None,                 # Controls where in foreground/background data/lines/bands appear.

    # Options for the legend.
    # 'upper right'   : 1
    # 'upper left'    : 2
    # 'lower left'    : 3
    # 'lower right'   : 4
    # 'right'         : 5
    # 'center left'   : 6
    # 'center right'  : 7
    # 'lower center'  : 8
    # 'upper center'  : 9
    # 'center'        : 10
    # 'best'          : Tries its best to automatically find a place for the legend.
    'legendpos': 'best',
    'bbox_to_anchor': None,      # Manual position of the legend.
    'legend_ncol': 1,            # Number of columns in the legend.
    'legend_col_spacing': None,  # Spacing between columns in the legend.
    'handletextpad': 0.2,        # Spacing between symbol and text in legend.
    'legend_title': None,        # Title of the legend.
    'alpha_legend': 0,           # Transperancy for the legend

    # Options for plotting files and functions.
    'xcol': None,        # Column for the xdata when plotting a file.
    'ycol': None,        # Column for the ydata when plotting a file.
    'yecol': None,       # Column for the errors in y-direction when plotting a file.
    'xecol': None,       # Column for the errors in x-direction when plotting a file.
    'style': "dots",     # Style when plotting a file.
    'npoints': 1000,     # Number of points for function plotting

    # Adjust aspects of the plot's axes
    'xtick_freq': None,
    'ytick_freq': None,
    'xtick_every_n' : 2,   # If xtick_freq or ytick_freq is not None, label every nth tick.
    'ytick_every_n' : 2,
    'xmin': None,          # Does not directly change x-range
    'xmax': None,
    'ymin': None,          # Does not directly change y-range
    'ymax': None,
}


# Used for later checks that the user did not pass a wrong parameter by mistake.
allowed_params = set(default_params.keys())
allowed_params = allowed_params | {'linestyle','show_leg','ha','va'}


# ---------------------------------------------------------------------------------------------- SOME EXTERNAL FUNCTIONS


def getColorGradient(NUM_COLORS=None):
    """ Return a preceptually uniform set of NUM_COLORS colors. """
    if NUM_COLORS is None:
        logger.TBError("Give me a number of colors.")
    cm = plt.get_cmap('viridis')
    gradColors=[]
    for i in range(NUM_COLORS):
        color = cm(1.*i/NUM_COLORS)
        gradColors.append(color)
    return gradColors


def latexify(bold=False):
    """ Allows use of LaTeX symbols in plots. The physics package is included, allowing use of
        convenient functions like ev. """
    logger.debug("Using latexify can be slow. If you need to speed up your plotting, suppress it.")
    if bold:
        plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{physics}\boldmath"
    else:
        plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{physics}"
    plt.rcParams['text.usetex'] = True


# ---------------------------------------------------------------------------------------------- SOME INTERNAL FUNCTIONS


def initializePlt(size,xmin,xmax,ymin,ymax):
    global initialize
    if initialize:
        initialize = False
        plt.rcParams['figure.autolayout'] = True
        plt.rcParams['axes.titlesize'] = size
        plt.rcParams['savefig.bbox'] = 'standard'
        plt.rcParams['ytick.labelsize'] = size
        plt.rcParams['font.size'] = size
        plt.rcParams['axes.labelsize'] = size
        plt.rcParams['legend.fontsize'] = size
        plt.rcParams['xtick.labelsize'] = size
        plt.rcParams['font.weight'] = default_params['font_weight']
        set_xrange(xmin,xmax)
        set_yrange(ymin,ymax)


def save_func(func, filename, args=(), func_err=None, args_err=(), grad = None, func_sup_numpy = False, **params):
    fill_param_dict(params)
    xmin = params['xmin']
    xmax = params['xmax']

    if params['expand']:
        wrap_func = lambda x, *wrap_args: func(x, *wrap_args)
        wrap_func_err = lambda x, *wrap_args_err: func_err(x, *wrap_args_err)
        wrap_grad = lambda x, *wrap_args: grad(x, *wrap_args)
    else:
        wrap_func = lambda x, *wrap_args: func(x, wrap_args)
        wrap_func_err = lambda x, *wrap_args_err: func_err(x, wrap_args_err)
        wrap_grad = lambda x, *wrap_args: grad(x, wrap_args)

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

        # Optional arguments that are constant and, therefore, not part of the error propagation
        tmp_opt = tuple(args)[len(args_err):]

        if func_sup_numpy:
            ydata_err = error_prop_func(xdata, wrap_func, tmp_args, args_err, grad = wrap_grad, args = tmp_opt)
        else:
            ydata_err = np.array([error_prop_func(x, wrap_func, tmp_args, args_err, grad = wrap_grad,
                                                  args = tmp_opt) for x in xdata])

        with open(filename, "w") as fout:
            for i in range(len(xdata)):
                print(xdata[i], ydata[i], ydata_err[i], file = fout)

    else:
        with open(filename, "w") as fout:
            for i in range(len(xdata)):
                print(xdata[i], ydata[i], file = fout)


def remove_points(data, *args, minval = None, maxval = None):
    if minval is None:
        minval = -np.inf
    if maxval is None:
        maxval = np.inf
    ind = (data>=minval) & (data<=maxval)
    ret = [data[ind]]
    for i in args:
        try:
            if i is not None:
                ret.append(i[ind])
            else:
                ret.append(None)
        except (IndexError, TypeError):
            ret.append(i)
    return ret


def clear_legend_labels():
    global legend_labels
    legend_labels = []
    global legend_handles
    legend_handles = []


def get_legend_handles():
    return legend_handles, legend_labels


def set_markers(marker_set=None):
    if marker_set is None:
        marker_set = markers_1
    global markers
    markers = itertools.cycle(marker_set)


def check_numpy(*data):
    data=list(data)
    for i in range(len(data)):
        if isinstance(data[i], (list, tuple)):
            data[i] = np.array(data[i])
    return tuple(data)


def getAxObject(params):
    if params['ax']==plt:
        return plt.gca()
    else:
        return params['ax']


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


def set_default_param(**kwargs):
    for key, val in kwargs.items():
        default_params[key] = val


def fill_param_dict(params):
    """Collection of default parameters for plotting routines. If a key does not exist in params, it is defined with
    a default value.
    
        Parameters
        ----------
        params: dictionary
            Dictionary with all parameters that are already set by the user

        Returns
        -------
            Dictionary filled with all default parameters
    """
    for key in params:
        if not key in allowed_params:
            logger.warn("Encountered unexpected plotting parameter",key)

    if 'show_leg' not in params:
        # When filling params for the first time, check if we show the legend. Triggered by one of the following keys
        for key in ('legend_title', 'legendpos', 'legend_ncol', 'legendpos_col_spacing', 'label', 'alpha_legend'):
            if key in params:
                params['show_leg'] = True
    if 'show_leg' not in params:
        params['show_leg'] = False

    for key, val in default_params.items():
        params.setdefault(key,val)


def add_optional(params):
    """Optional parameter that are not defined in fill_param_dict are collected by this function.

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


def set_params(**params):
    """Set additional parameters to the plot. For example set a title or label.
        Parameters
        ----------

        **params :
            Additional parameters that can be set.
    """

    logger.debug("Note that many of the plotting functions call set_params, which can reset what you pass as argument.")
    logger.debug("If you are having trouble passing options to set_params, try calling it at the end of your script.")

    fill_param_dict(params)

    initializePlt(params['font_size'],params['xmin'],params['xmax'],params['ymin'],params['ymax'])

    ax  = getAxObject(params)
    zod = params['zod']

    if zod is None:
         zod = globals()['zod']

    if params['xlabel'] is not None:
        if params['labelsintoplot'] or params['xlabelpos'] is not None:
            if params['xlabelpos'] is None:
                params['xlabelpos'] = (0.95,0.027)
            ax.annotate(params['xlabel'], xy=params['xlabelpos'], xycoords='axes fraction', color='black',
                        fontsize=params['font_size'], fontweight=params['font_weight'], ha='right', va='bottom',
                        bbox=dict(linewidth=0, facecolor='white', alpha=params['alpha_label']), zorder=zod)
        else:
            ax.xlabel(params['xlabel'])
            ax.xaxis.get_label().set_fontsize(params['font_size'])

    if params['ylabel'] is not None:
        if params['labelsintoplot'] or params['ylabelpos'] is not None:
            if params['ylabelpos'] is None:
                params['ylabelpos'] = (0.025,0.972)
            ax.annotate(params['ylabel'], xy=params['ylabelpos'], xycoords='axes fraction', color='black',
                        fontsize=params['font_size'], ha='left', va='top', fontweight=params['font_weight'],
                        bbox=dict(linewidth=0, facecolor='white', alpha=params['alpha_label']), zorder=zod)
        else:
            ax.ylabel(params['ylabel'])
            ax.yaxis.get_label().set_fontsize(params['font_size'])

    if params['title'] is not None:
        plt.title(params['title'])

    if params['ticksintoplot'] is not None:
        ax.tick_params(direction='in')

    if params['surroundWithTicks']:
        ax.tick_params(top=True,right=True)

    if params['show_leg']:
        leg = ax.legend(legend_handles, legend_labels, numpoints=1, bbox_to_anchor = params['bbox_to_anchor'],
                        title=params['legend_title'], loc=params['legendpos'], ncol=params['legend_ncol'],
                        columnspacing=params['legend_col_spacing'],handletextpad = params['handletextpad'])
        leg.get_frame().set_alpha(params['alpha_legend'])
        leg.set_zorder(zod)

    if params['xtick_freq'] is not None:
        start, end = ax.get_xlim()
        if params['xmin'] is not None:
            start = params['xmin']
        if params['xmax'] is not None:
            end = params['xmax']
        ax.xaxis.set_ticks(np.arange(start, end, params['xtick_freq']))
        for n, label in enumerate(ax.xaxis.get_ticklabels()):
            if n % params['xtick_every_n'] != 0:
                label.set_visible(False)

    if params['ytick_freq'] is not None:
        start, end = ax.get_ylim()
        if params['ymin'] is not None:
            start = params['ymin']
        if params['ymax'] is not None:
            end = params['ymax']
        ax.yaxis.set_ticks(np.arange(start, end, params['ytick_freq']))
        for n, label in enumerate(ax.yaxis.get_ticklabels()):
            if n % params['ytick_every_n'] != 0:
                label.set_visible(False)


# ------------------------------------------------------------------------------------------------ MAIN PLOTTING METHODS


def plot_file(filename, xcol=1, ycol=2, yecol=None, xecol=None, func = None, func_args = (), **params):
    fill_param_dict(params)
    data = np.loadtxt(filename, dtype = np.str).transpose()
    if xcol is not None:
        xdata = np.array(data[xcol - 1], dtype = float)
    else:
        xdata = np.arange(len(data))
    ydata = np.array(data[ycol - 1], dtype = float)
    yedata = None
    xedata = None

    if yecol is not None:
        yedata = np.array(data[yecol - 1], dtype = float)
    else:
        if params['style'] == "fill":
            logger.TBError("Need error column for filled plotting")
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

    logger.TBError("No such style: " + params['style'])


def plot_dots(xdata, ydata, yedata = None, xedata = None, **params):
    xdata, ydata, yedata, xedata = check_numpy(xdata, ydata, yedata, xedata)
    fill_param_dict(params)
    optional = add_optional(params)

    xdata, ydata, yedata, xedata = remove_points(xdata, ydata, yedata, xedata, minval=params['xmin'], maxval=params['xmax'])

    ax  = getAxObject(params)

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
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], alpha=params['alpha_dots'],  color=params['color'],
                           zorder=zod, markersize=params['markersize'], capsize=params['capsize'],
                           elinewidth=params['elinewidth'], markeredgewidth=params['elinewidth'],
                           markerfacecolor = params['point_fill_color'], **optional)
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], alpha=params['alpha_dots'], zorder=zod,
                           markersize=params['markersize'], capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], markerfacecolor = params['point_fill_color'],
                           **optional)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append(ebar)

    globals()['zod'] += 1
    set_params(**params)
    
    return ebar


def plot_bar(xdata, ydata, width=None, align='edge', alpha=1.0, edgecolor='#666677',linewidth=0.2, **params):
    if width is None:
        width = xdata[1] - xdata[0]

    if alpha is None:
        alpha=params['alpha']

    xdata, ydata = check_numpy(xdata, ydata)
    fill_param_dict(params)
    optional = add_optional(params)
    ax  = getAxObject(params)

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
    xdata, ydata, yedata, xedata = remove_points(xdata, ydata, yedata, xedata, minval=params['xmin'], maxval=params['xmax'])
    xscale=params['xscale']
    yscale=params['yscale']
    ax  = getAxObject(params)

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
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], color=params['color'], zorder=zod,
                           markersize=params['markersize'], capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], alpha=params['alpha_dots'],
                           markerfacecolor = params['point_fill_color'])
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], zorder=zod, markersize=params['markersize'],
                           capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], alpha=params['alpha_dots'],
                           markerfacecolor = params['point_fill_color'])

    globals()['zod'] += 1

    col = ebar[0].get_color()

    line = ax.errorbar(xdata*xscale, ydata*yscale, color = col, linewidth=params['linewidth'], zorder = zod,
                       alpha = params["alpha_lines"], **optional)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append((line, ebar))

    set_params(**params)
    return ebar


def plot_func(func, args=(), func_err=None, args_err=(), grad = None, func_sup_numpy = False, swapXY=False, **params):
    """ To plot an error band with an explicit error function, use func_err. args_err are all parameters for func_err.
        To use a numerical derivative, just pass the errors of args to args_err. The option swapXY allows for error
        bars in the x-direction instead of the y-direction. """

    fill_param_dict(params)
    params['marker'] = None
    xmin = params['xmin']
    xmax = params['xmax']

    if params['expand']:
        wrap_func = lambda x, *wrap_args: func(x, *wrap_args)
        wrap_func_err = lambda x, *wrap_args_err: func_err(x, *wrap_args_err)
        wrap_grad = lambda x, *wrap_args: grad(x, *wrap_args)
    else:
        wrap_func = lambda x, *wrap_args: func(x, wrap_args)
        wrap_func_err = lambda x, *wrap_args_err: func_err(x, wrap_args_err)
        wrap_grad = lambda x, *wrap_args: grad(x, wrap_args)

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

        if swapXY:
            return plot_fill(xdata, ydata, yedata=None, xedata=ydata_err, **params)
        else:
            return plot_fill(xdata, ydata, ydata_err, **params)

    elif len(args_err) > 0:
        if grad is None:
            logger.warn("Used numerical derivative!")
            wrap_grad = None

        # Arguments that are part of the error propagation
        tmp_args = tuple(args)[0:len(args_err)]

        # Optional arguments that are constant and, therefore, not part of the error propagation
        tmp_opt = tuple(args)[len(args_err):]

        if func_sup_numpy:
            ydata_err = error_prop_func(xdata, wrap_func, tmp_args, args_err, grad = wrap_grad, args = tmp_opt)
        else:
            ydata_err = np.array([error_prop_func(x, wrap_func, tmp_args, args_err, grad = wrap_grad,
                                                  args = tmp_opt) for x in xdata])

        if swapXY:
            return plot_fill(xdata, ydata, yedata=None, xedata=ydata_err, **params)
        else:
            return plot_fill(xdata, ydata, ydata_err, **params)

    else:
        if swapXY:
            return plot_lines(ydata, xdata, yedata=None, xedata=None, **params)
        else:
            return plot_lines(xdata, ydata, yedata=None, xedata=None, **params)


def plot_fill(xdata, ydata, yedata, xedata=None, pattern=None, **params):

    if (yedata is None) and (xedata is None):
        logger.TBError("Please pass plot_fill some error bars.")
    if (yedata is not None) and (xedata is not None):
        logger.TBError("Please pass plot_fill either x-error or y-error, not both.")

    fill_param_dict(params)
    optional = add_optional(params)

    if xedata is None:
        xdata, ydata, yedata = check_numpy(xdata, ydata, yedata)
        xdata, ydata, yedata = remove_points(xdata, ydata, yedata, minval=params['xmin'], maxval=params['xmax'])
    else:
        xdata, ydata, xedata = check_numpy(xdata, ydata, xedata)
        ydata, xdata, xedata = remove_points(ydata, xdata, xedata, minval = params['ymin'], maxval = params['ymax'])

    xscale = params['xscale']
    yscale = params['yscale']
    ax  = getAxObject(params)
    zod    = params['zod']
    if zod is None:
         zod = globals()['zod']

    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, linewidth=params['linewidth'], color=params['color'], zorder=zod,
                           alpha = params['alpha_lines'], **optional)
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, linewidth=params['linewidth'], zorder=zod,
                           alpha = params['alpha_lines'], **optional)

    globals()['zod'] += 1

    col = ebar[0].get_color()
    if xedata is None:
        pl = ax.fill_between(xdata*xscale, (np.asarray(ydata*yscale) - np.asarray(yedata*yscale)),
                             (np.asarray(ydata*yscale) + np.asarray(yedata*yscale)), facecolor=col, alpha=params['alpha'],
                             linewidth=0, zorder=1, hatch=pattern, edgecolor=col)
    else:
        pl = ax.fill_betweenx(ydata * yscale, (np.asarray(xdata * xscale) - np.asarray(xedata * xscale)),
                             (np.asarray(xdata * xscale) + np.asarray(xedata * xscale)), facecolor=col,
                             alpha=params['alpha'],linewidth=0, zorder=1, hatch=pattern, edgecolor=col)

    if params['label'] is not None:
        legend_labels.append(params['label'])
        legend_handles.append((ebar, pl))

    set_params(**params)
    return ebar,pl


#
# TODO: make this guy work vertically
#
def plot_band(xdata, low_lim, up_lim, center = None, **params):
    fill_param_dict(params)
    optional = add_optional(params)
    if center is not None:
        xdata, low_lim, up_lim, center = remove_points(xdata, low_lim, up_lim, center, minval = params['xmin'], maxval = params['xmax'])
    else:
        xdata, low_lim, up_lim = remove_points(xdata, low_lim, up_lim, minval = params['xmin'], maxval = params['xmax'])

    xscale=params['xscale']
    yscale=params['yscale']
    ax  = getAxObject(params)

    zod = params['zod']
    if zod is None:
         zod = globals()['zod']

    globals()['zod'] += 1

    if params['color'] is None:
        pl = ax.fill_between(xdata*xscale, yscale*low_lim, yscale*up_lim, alpha=params['alpha'], linewidth=0, zorder=1)
    else:
        pl = ax.fill_between(xdata*xscale, yscale*low_lim, yscale*up_lim, facecolor=params['color'],
                             alpha=params['alpha'], linewidth=0, zorder=1)

    col = cl.rgb2hex(pl.get_facecolor()[0])

    if params['alpha_lines'] != 0:
        ax.errorbar(xdata*xscale, yscale*low_lim, color = col, linewidth=params['linewidth'], zorder = zod,
                    alpha = params["alpha_fill_edge"])
        ax.errorbar(xdata*xscale, yscale*up_lim, color = col, linewidth=params['linewidth'], zorder = zod,
                    alpha = params["alpha_fill_edge"])

    ebar = None
    if center is not None:
        ebar = ax.errorbar(xdata*xscale, center*yscale, linewidth=params['linewidth'], color=col, zorder=zod,
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


def plot_cov(cov, filename = None, title=None, ignore_first = 0, norm = True, xrange = None, yrange = None,
             xmin = None, xmax = None, ymin = None, ymax = None, xlabel = "$n_{\\tau/\\sigma}$",
             ylabel = "$n_{\\tau/\\sigma}$"):
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


def plot_eig(cov, filename, title=None):
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