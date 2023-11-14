#
# plotting.py
#
# H. Sandmeyer, D. Clarke
#
# Collection of convenience tools for plotting using matplotlib.
#

import itertools, warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import matplotlib.ticker as ticker
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkEqualLengths, checkType
from latqcdtools.base.utilities import isHigherDimensional, toNumpy
from latqcdtools.base.readWrite import readTable
warnings.filterwarnings("ignore", category=UserWarning)


ZOD        = 1      # Orders different layers 
INITIALIZE = True   # A global flag to ensure we only initialize once
LEGEND     = False  # A global flag to apply legend attributes when we have one
FOREGROUND = 99999  # zorder = FOREGROUND to put something in plot on top of everything else 
BACKGROUND = 0      # zorder = BACKGROUND to put something in plot behind everything else 


colors = ['#d32d11', '#0081bf', '#e5af11', '#7c966d', '#7570b3', '#ff934f', '#666666', '#D186B3']


markers_1 = ['o', 'v', 'D', 's', 'p', '^', 'h', 'd', 'x', '+', '*']
markers   = itertools.cycle(markers_1)


# Handles are what graphics you're going to use to represent a particular legend element. Usually
# these are taken from some feature of the plotted data. The labels are the names they will have.
# We implement these as dictionaries to allow the user to have different legend objects in
# for instance multiple subplots.
legend_handles = { plt : [] }
legend_labels  = { plt : [] }


default_params = {

    'ax': plt,  # Axis object that is used for the plots. Default is matplotlib.pyplot.

    # Basic options affecting most plots.
    'xlabel': None,
    'alpha_xlabel' : 1,          # Transparency for x-label
    'ylabel': None,
    'alpha_ylabel' : 1,          # Transparency for y-label
    'xmin': None,
    'xmax': None,
    'ymin': None,
    'ymax': None,
    'title': None,
    'label': None,               # What are the data called? (Will appear in legend.)
    'color': None,               # Color for your data. (By default each new set automatically gets different color.)
    'marker': "iter",            # Symbol used for plotting data. (Set to 'None' if you don't want any.)
    'markersize': 8,             # Size of the symbols.
    'font_size': 16,             # Default font size for text.
    'font_weight': 'normal',     # Default style of font ('normal', 'bold', 'heavy', 'light')
    'alpha': 0.5,                # General transparency for data.
    'ticksintoplot': True,       # Put ticks into plotting area.
    'surroundWithTicks': True,   # Put ticks also on top and right.
    'labelsintoplot': True,      # Put xlabel and ylabel into plotting area.
    'xlabelpos': None,           # If labelsintplot=True, shift the position (x,y) of the x-label, expressed as percent.
    'ylabelpos': None,
    'alpha_dots': None,          # Transperancy for different dots
    'alpha_lines': 1,            # Transperancy for different lines
    'alpha_fill_edge': 0,        # Transperancy for edges of error bands
    'linewidth': 1,              # Linewidth of line plots
    'capsize': 1.5,              # Length of caps af error bars
    'elinewidth': 1.0,           # Linewidth of the error bars of caps af error bars
    'point_fill_color': "None",  # Fill color of points. Set to None (not as string) to have filled symbols
    'ZOD': None,                 # Controls where in foreground/background data/lines/bands appear.

    # Options for the legend.
    # legendpos:
    #   2     9     1
    #   6     10    7
    #   3     8     4
    'legendpos': 'best',
    'bbox_to_anchor': None,      # Manual position of the legend. The very bottom-left is (0,0), and the very 
                                 #   top-right is (1,1). If you set this, legendpos appears to get ignored.
    'legend_ncol': 1,            # Number of columns in the legend.
    'legend_col_spacing': None,  # Spacing between columns in the legend.
    'handletextpad': 0.2,        # Spacing between symbol and text in legend.
    'legend_title': None,        # Title of the legend.
    'alpha_legend': 1,           # Transperancy for the legend.

    # Adjust special aspects of the plot's axes
    'xtick_freq': None,
    'ytick_freq': None,
    'xtick_every_n': 2,    # If xtick_freq or ytick_freq is not None, label every nth tick.
    'ytick_every_n': 2,
    'xtick_format' : None, # Format the y-ticks, e.g. if you want to specify the number of decimals. 
    'ytick_format' : None,
    'xscale': 1.0,         # Scale data in xdata by this factor.
    'yscale': 1.0,
    'xlogscale': False,    # Should we use a log scale for the x-axis?
    'ylogscale': False,
}


# Used for later checks that the user did not pass a wrong parameter by mistake.
allowed_params = set(default_params.keys())
allowed_params = allowed_params | {'linestyle','ha','va'}


# ---------------------------------------------------------------------------------------------- SOME EXTERNAL FUNCTIONS


def latexify(bold=False):
    """ Allows use of LaTeX symbols in plots. The physics package is included, allowing use of
        convenient functions like ev. """
    logger.info("Using LaTeX to make pretty plots.")
    if bold:
        plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{physics}\boldmath"
    else:
        plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{physics}"
    plt.rcParams['text.usetex'] = True


def clearPlot():
    """ Clears plot object, legend handles, and zorder. Useful if you want to do multiple plots in the same script. """
    logger.debug('Reset plot defaults.')
    global INITIALIZE, ZOD, LEGEND, legend_labels, legend_handles
    ZOD        = 1
    INITIALIZE = True
    LEGEND     = False
    legend_handles = { plt : [] }
    legend_labels  = { plt : [] }
    plt.clf()


def getColorGradient(NUM_COLORS,map='viridis'):
    """ Generate perceptually uniform set of colors. Useful when you need more than 8 colors.

    Args:
        NUM_COLORS (int): number of colors you need 
        map (str, optional): use custom colormap. Defaults to 'viridis'.

    Returns:
        list: colors 
    """
    checkType(NUM_COLORS,int)
    checkType(map,str)
    cm = plt.get_cmap(map)
    gradColors=[]
    for i in range(NUM_COLORS):
        color = cm(1.*i/NUM_COLORS)
        gradColors.append(color)
    return gradColors


def set_xrange(xmin=None, xmax=None, ax=plt):
    _set_xmin(ax,xmin)
    _set_xmax(ax,xmax)


def set_yrange(ymin=None, ymax=None, ax=plt):
    _set_ymin(ax,ymin)
    _set_ymax(ax,ymax)


def set_default_param(**kwargs):
    """ Lets the user adjust the default parameter settings. """
    for key, val in kwargs.items():
        default_params[key] = val


# ---------------------------------------------------------------------------------------------- SOME INTERNAL FUNCTIONS


def _initializePlt(size,xmin,xmax,ymin,ymax):
    """ Set up inital plot parameters, like its size. I tried to introduce a global variable INITIALIZE that checks
        so that this only gets called once per plot script. """
    global INITIALIZE
    if INITIALIZE:
        logger.debug("Plot initializer called!")
        logger.debug("Many of the plotting functions call set_params, which can reset what you pass as argument.")
        logger.debug("If you have trouble passing options to set_params, try calling it at the end of your script.")
        INITIALIZE = False
        plt.rcParams['figure.autolayout'] = True
        plt.rcParams['axes.titlesize'] = size
        plt.rcParams['savefig.bbox'] = 'standard'
        plt.rcParams['font.size'] = size
        plt.rcParams['ytick.labelsize'] = size
        plt.rcParams['xtick.labelsize'] = size
        plt.rcParams['axes.labelsize'] = size
        plt.rcParams['legend.fontsize'] = size
        plt.rcParams['font.weight'] = default_params['font_weight']
        set_xrange(xmin,xmax)
        set_yrange(ymin,ymax)


def _remove_points(data, *args, minval=None, maxval=None):
    """ Assuming *args are a tuple of arrays of the same length as data, cut data and
    args so that data fall between minval and maxval, and keep only the corresponding
    elements of args. 

    Args:
        data (array-like)
        minval (float, optional): minimum allowed data. Defaults to -np.inf.
        maxval (float, optional): maximum allowed data. Defaults to np.inf.

    Returns:
        tuple: data, *args trimmed according to minval and maxval
    """
    checkType(data,"array")
    data = np.array(data)
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


def _update_labels(ax,label):
    checkType(label,str)
    global legend_labels
    if ax in legend_labels:
        legend_labels[ax].append(label)
    else:
        legend_labels[ax] = [label]


def _update_handles(ax,handle):
    global legend_handles
    if ax in legend_handles:
        legend_handles[ax].append(handle)
    else:
        legend_handles[ax] = [handle]


def _getAxObject(params):
    if params['ax']==plt:
        return plt.gca()
    else:
        return params['ax']


def _set_xmin(ax,x_min=None):
    if x_min is not None:
        x1, x2 = ax.get_xlim()
        ax.set_xlim([x_min,x2])


def _set_xmax(ax,x_max=None):
    if x_max is not None:
        x1, x2 = ax.get_xlim()
        ax.set_xlim([x1,x_max])


def _set_ymin(ax,y_min=None):
    if y_min is not None:
        y1, y2 = ax.get_ylim()
        ax.set_ylim([y_min,y2])


def _set_ymax(ax,y_max=None):
    if y_max is not None:
        y1, y2 = ax.get_ylim()
        ax.set_ylim([y1,y_max])


def fill_param_dict(params):
    """ Collection of default parameters for plotting routines. If a key does not exist in params, it is defined with
    a default value.
    
        Parameters
        ----------
        params: dictionary
            Dictionary with all parameters that are already set by the user

    """
    global LEGEND
    for key in params:
        if not key in allowed_params:
            logger.warn("Encountered unexpected plotting parameter",key)

    if not LEGEND:
        # When filling params, check if we show the legend. This is triggered by one of these keys
        # being different from the default value. 
        for key in ('legend_title', 'legendpos', 'legend_ncol', 'legendpos_col_spacing', 'label', 'alpha_legend'):
            if key in params:
                if params[key] != default_params[key]:
                    logger.debug('Found legend trigger',key)
                    LEGEND = True

    for key, val in default_params.items():
        params.setdefault(key,val)


def _add_optional(params):
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
    ret = {}
    for key in params:
        if not key in reference:
            ret[key] = params[key]
    return ret


def set_params(**params):
    """ Set additional parameters to the plot. For example set a title or label.

        Args:
            **params: Additional parameters that can be set.
    """
    global LEGEND

    fill_param_dict(params)

    _initializePlt(params['font_size'],params['xmin'],params['xmax'],params['ymin'],params['ymax'])

    ax  = _getAxObject(params)
    ZOD = params['ZOD']

    if ZOD is None:
        ZOD = globals()['ZOD']

    if params['xlabelpos'] is not None:
        if params['xlabel'] is None:
            logger.warn('Set xlabelpos with no xlabel.')
    if params['ylabelpos'] is not None:
        if params['ylabel'] is None:
            logger.warn('Set ylabelpos with no ylabel.')

    if params['xlabel'] is not None:
        checkType(params['xlabel'],str)
        if params['labelsintoplot'] or params['xlabelpos'] is not None:
            if params['xlabelpos'] is None:
                params['xlabelpos'] = (0.95,0.027)
            ax.annotate(params['xlabel'], xy=params['xlabelpos'], xycoords='axes fraction', color='black',
                        fontsize=params['font_size'], fontweight=params['font_weight'], ha='right', va='bottom',
                        bbox=dict(linewidth=0, facecolor='white', edgecolor=None, alpha=params['alpha_xlabel']), 
                        zorder=FOREGROUND)
        else:
            ax.set_xlabel(params['xlabel'])
            ax.xaxis.get_label().set_fontsize(params['font_size'])

    if params['ylabel'] is not None:
        checkType(params['ylabel'],str)
        if params['labelsintoplot'] or params['ylabelpos'] is not None:
            if params['ylabelpos'] is None:
                params['ylabelpos'] = (0.025,0.972)
            ax.annotate(params['ylabel'], xy=params['ylabelpos'], xycoords='axes fraction', color='black',
                        fontsize=params['font_size'], ha='left', va='top', fontweight=params['font_weight'],
                        bbox=dict(linewidth=0, facecolor='white', edgecolor=None, alpha=params['alpha_ylabel']), 
                        zorder=FOREGROUND)
        else:
            ax.set_ylabel(params['ylabel'])
            ax.yaxis.get_label().set_fontsize(params['font_size'])

    if params['title'] is not None:
        checkType(params['title'],str)
        plt.title(params['title'])

    if params['xlogscale'] and params['xtick_freq'] is not None:
        logger.warn("xtick_freq assumes no log scale.")

    if params['ylogscale'] and params['ytick_freq'] is not None:
        logger.warn("ytick_freq assumes no log scale.")

    if params['xlogscale']:
        ax.set_xscale('log')

    if params['ylogscale']:
        ax.set_yscale('log')

    if params['ticksintoplot'] is not None:
        ax.tick_params(which='both',direction='in')

    if params['surroundWithTicks']:
        ax.tick_params(top=True,right=True)

    set_xrange(params['xmin'],params['xmax'],ax)
    set_yrange(params['ymin'],params['ymax'],ax)

    if LEGEND:
        if not ax in legend_handles:
            logger.TBError('Legend for axis',ax,'was activated without any label.')
        leg = ax.legend(legend_handles[ax], legend_labels[ax], numpoints=1, bbox_to_anchor = params['bbox_to_anchor'],
                        title=params['legend_title'], loc=params['legendpos'], ncol=params['legend_ncol'],
                        columnspacing=params['legend_col_spacing'],handletextpad = params['handletextpad'])
        leg.get_frame().set_alpha(params['alpha_legend'])
        leg.set_zorder(FOREGROUND)
        plt.tight_layout()

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

    if params['xtick_format'] is not None:
        checkType(params['xtick_format'],str)
        if plt.rcParams['text.usetex']:
            x_format = ticker.StrMethodFormatter('$'+params['xtick_format']+'$')
        else:
            x_format = ticker.StrMethodFormatter(params['xtick_format'])
        ax.get_xaxis().set_major_formatter(x_format)

    if params['ytick_format'] is not None:
        checkType(params['ytick_format'],str)
        if plt.rcParams['text.usetex']:
            y_format = ticker.StrMethodFormatter('$'+params['ytick_format']+'$')
        else:
            y_format = ticker.StrMethodFormatter(params['ytick_format'])
        ax.get_yaxis().set_major_formatter(y_format)


# ------------------------------------------------------------------------------------------------ MAIN PLOTTING METHODS


def preliminary(x,y,text='PRELIMINARY',**kwargs):
    """ Generate a PRELIMINARY tag on the plot.

    Args:
        x (float): x-position of bottom-left corner (in units of x-axis) 
        y (float): y-position of bottom-left corner (in units of y-axis) 
        text (str, optional): Text indicating result is preliminary. Defaults to 'PRELIMINARY'.
    """
    checkType(text,str)
    if 'color' in kwargs: 
        color=kwargs['color']
    else:
        color='gray'
    plt.text(x, y, text, color=color)


def plot_file(filename, xcol=0, ycol=1, yecol=None, xecol=None, func = None, func_args = (), style='dots', **params):
    """ Plot data in file. You can set the style with the style argument. Columns indexed from 0.

    Args:
        filename (str): _description_
        xcol (int, optional): Which column is xdata. Defaults to 1.
        ycol (int, optional): Which column is ydata. Defaults to 2.
        yecol (int, optional): Which column has y error. Defaults to None.
        xecol (int, optional): Which column has x error. Defaults to None.
        func (function, optional): Apply this function to the data. Defaults to None.
        func_args (tuple, optional): Arguments to func. Defaults to ().
        style (str, optional): Choose from dots, lines, fill, and band. Defaults to 'dots'.
        **params: Additional parameters that can be set.
    """
    fill_param_dict(params)
    data = readTable(filename,dtype=str)
    checkType(xcol,int)
    checkType(ycol,int)
    xdata  = data[xcol].astype(float)
    ydata  = data[ycol].astype(float)
    yedata = None
    xedata = None

    if yecol is not None:
        yedata = data[yecol].astype(float)
    else:
        if style == "fill":
            logger.TBError("Need error column for filled plotting")
    if xecol is not None:
        xedata = data[xecol].astype(float)

    checkEqualLengths(xdata,ydata,xedata,yedata)

    if func is not None:
        if yedata is not None:
            if xedata is not None:
                xdata, ydata, yedata, xedata = func(xdata, ydata, yedata, xedata, *func_args)
            else:
                xdata, ydata, yedata = func(xdata, ydata, yedata, *func_args)
        else:
            xdata, ydata = func(xdata, ydata, *func_args)

    if style == "dots":
        return plot_dots(xdata, ydata, yedata=yedata, xedata=xedata, **params)
    elif style == "lines":
        return plot_lines(xdata, ydata, yedata=yedata, xedata=xedata, **params)
    elif style == "fill":
        return plot_fill(xdata, ydata, yedata=yedata, **params)
    elif style == "band":
        return plot_band(xdata, ydata, yedata, xedata, **params)
    else:
        logger.TBError("Unknown style",style)


def plot_dots(xdata, ydata, yedata = None, xedata = None, **params):
    """ Plot ydata vs xdata as dots. 

    Args:
        xdata (array-like)
        ydata (array-like)
        yedata (array-like, optional): y error. Defaults to None.
        xedata (array-like, optional): x error. Defaults to None.
        **params: Additional parameters that can be set.
    """
    checkEqualLengths(xdata,ydata,xedata,yedata)
    xdata, ydata, xedata, yedata = toNumpy(xdata, ydata, xedata, yedata)
    fill_param_dict(params)
    optional = _add_optional(params)

    xdata, ydata, yedata, xedata = _remove_points(xdata, ydata, yedata, xedata, minval=params['xmin'], maxval=params['xmax'])

    ax  = _getAxObject(params)

    xscale = params['xscale']
    yscale = params['yscale']
    if xedata is not None:
        xedata=np.copy(xedata*xscale)
    if yedata is not None:
        yedata=np.copy(yedata*yscale)

    ZOD = params['ZOD']
    if ZOD is None:
        ZOD = globals()['ZOD']

    marker = params['marker']
    if marker == "iter":
        marker = next(markers)
    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], alpha=params['alpha_dots'],  color=params['color'],
                           zorder=ZOD, markersize=params['markersize'], capsize=params['capsize'],
                           elinewidth=params['elinewidth'], markeredgewidth=params['elinewidth'],
                           markerfacecolor = params['point_fill_color'], **optional)
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], alpha=params['alpha_dots'], zorder=ZOD,
                           markersize=params['markersize'], capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], markerfacecolor = params['point_fill_color'],
                           **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,ebar)

    globals()['ZOD'] += 1
    set_params(**params)
    
    return ebar


def plot_bar(xdata, ydata, width=None, align='edge', alpha=1.0, edgecolor='#666677',linewidth=0.2, **params):
    """ Plot ydata vs xdata as bars.

    Args:
        xdata (array-like)
        ydata (array-like)
        width (_type_, optional): _description_. Defaults to None.
        align (str, optional): _description_. Defaults to 'edge'.
        alpha (float, optional): Transparency. Defaults to 1.0.
        edgecolor (str, optional): Color of bar edges. Defaults to '#666677'.
        linewidth (float, optional): Defaults to 0.2.
        **params: Additional parameters that can be set.
    """
    checkEqualLengths(xdata,ydata)
    xdata, ydata = toNumpy(xdata, ydata)
    fill_param_dict(params)
    optional = _add_optional(params)
    ax  = _getAxObject(params)

    if width is None:
        width = xdata[1] - xdata[0]

    if alpha is None:
        alpha=params['alpha']

    ZOD = params['ZOD']
    if ZOD is None:
        ZOD = globals()['ZOD']

    if params['color'] is not None:
        bar = ax.bar(xdata, ydata, color=params['color'], zorder=ZOD, width=width, align=align, edgecolor=edgecolor,
                     linewidth=linewidth, alpha=alpha, **optional)
    else:
        bar = ax.bar(xdata, ydata, zorder=ZOD, width=width, align=align, edgecolor=edgecolor,
                     linewidth=linewidth, alpha=alpha, **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,bar)

    globals()['ZOD'] += 1
    set_params(**params)

    return bar


def plot_hist(data, bins = None, density=False, label=None, **params):
    """ Create a histogram of the array data. If you would like to plot multiple data sets in the same histogram,
    simply pass as a list or tuple of arrays of data, like data = [list1, list2, ...].

    Args:
        data (array-like)
        bins (int, optional): Number of bins. Defaults to None, which sets the number of bins automatically.
        **params: Additional parameters that can be set.
    """
    fill_param_dict(params)
    ax = _getAxObject(params)
    if bins is None:
        bins = 'auto'
    if isHigherDimensional(data):
        ax.hist(data, bins=bins, density=density,label=label)
    else:
        ax.hist(data, bins=bins, density=density,label=label,color=params['color'])
    if density:
        ax.set_yticklabels([])
    if label is not None:
        ax.legend()
    set_params(**params)


def plot_lines(xdata, ydata, yedata=None, xedata=None, **params):
    """ Plot ydata vs xdata using lines.

    Args:
        xdata (array-like)
        ydata (array-like)
        yedata (array-like, optional): y error. Defaults to None.
        xedata (array-like, optional): x error. Defaults to None.
        **params: Additional parameters that can be set.
    """
    checkEqualLengths(xdata,ydata,yedata,xedata)
    xdata, ydata, xedata, yedata = toNumpy(xdata, ydata, xedata, yedata)
    fill_param_dict(params)
    optional = _add_optional(params)
    xdata, ydata, yedata, xedata = _remove_points(xdata, ydata, yedata, xedata, minval=params['xmin'], maxval=params['xmax'])
    xscale=params['xscale']
    yscale=params['yscale']
    ax  = _getAxObject(params)

    if xedata is not None:
        xedata=np.copy(xedata*xscale)
    if yedata is not None:
        yedata=np.copy(yedata*yscale)

    ZOD = params['ZOD']
    if ZOD is None:
        ZOD = globals()['ZOD']

    marker = params['marker']
    if marker == "iter":
        marker = next(markers)

    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], color=params['color'], zorder=ZOD,
                           markersize=params['markersize'], capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], alpha=params['alpha_dots'],
                           markerfacecolor = params['point_fill_color'])
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, yerr=yedata, xerr=xedata, marker=marker, linestyle='None',
                           linewidth=params['linewidth'], zorder=ZOD, markersize=params['markersize'],
                           capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], alpha=params['alpha_dots'],
                           markerfacecolor = params['point_fill_color'])

    globals()['ZOD'] += 1

    col = ebar[0].get_color()

    line = ax.errorbar(xdata*xscale, ydata*yscale, color = col, linewidth=params['linewidth'], zorder = ZOD,
                       alpha = params["alpha_lines"], **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,(line,ebar))

    set_params(**params)
    return ebar


def plot_fill(xdata, ydata, yedata, xedata=None, pattern=None, **params):
    """ Plot a filled region within ydata +/- yedata. Can set xedata along with yedata=None for vertical bands.

    Args:
        xdata (array-like)
        ydata (array-like)
        yedata (array-like): y error 
        xedata (array-like, optional): x error. Defaults to None.
        pattern (_type_, optional): _description_. Defaults to None.
        **params: Additional parameters that can be set.
    """
    if (yedata is None) and (xedata is None):
        logger.TBError("Please pass some error bars.")
    if (yedata is not None) and (xedata is not None):
        logger.TBError("Please pass either x-error or y-error, not both.")
    xdata, ydata, xedata, yedata = toNumpy(xdata, ydata, xedata, yedata)
    checkEqualLengths(xdata,ydata,xedata,yedata)

    fill_param_dict(params)
    optional = _add_optional(params)

    if xedata is None:
        xdata, ydata, yedata = _remove_points(xdata, ydata, yedata, minval=params['xmin'], maxval=params['xmax'])
    else:
        ydata, xdata, xedata = _remove_points(ydata, xdata, xedata, minval=params['ymin'], maxval=params['ymax'])

    xscale = params['xscale']
    yscale = params['yscale']
    ax     = _getAxObject(params)
    ZOD    = params['ZOD']
    if ZOD is None:
        ZOD = globals()['ZOD']

    if params['color'] is not None:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, linewidth=params['linewidth'], color=params['color'], zorder=ZOD,
                           alpha = params['alpha_lines'], **optional)
    else:
        ebar = ax.errorbar(xdata*xscale, ydata*yscale, linewidth=params['linewidth'], zorder=ZOD,
                           alpha = params['alpha_lines'], **optional)

    globals()['ZOD'] += 1

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
        _update_labels(ax,params['label'])
        _update_handles(ax,(ebar,pl))

    set_params(**params)
    return ebar,pl


def plot_band(xdata, low_lim, up_lim, center = None, **params):
    """ Plot a horizontal band.

    Args:
        xdata (array-like)
        low_lim (float): _description_
        up_lim (float): _description_
        center (_type_, optional): _description_. Defaults to None.
        **params: Additional parameters that can be set.
    """
    xdata, low_lim, up_lim, center = toNumpy(xdata, low_lim, up_lim, center)
    checkEqualLengths(xdata, low_lim, up_lim, center)
    fill_param_dict(params)
    optional = _add_optional(params)
    if center is not None:
        xdata, low_lim, up_lim, center = _remove_points(xdata, low_lim, up_lim, center, minval=params['xmin'], maxval=params['xmax'])
    else:
        xdata, low_lim, up_lim = _remove_points(xdata, low_lim, up_lim, minval=params['xmin'], maxval=params['xmax'])

    xscale=params['xscale']
    yscale=params['yscale']
    ax  = _getAxObject(params)

    ZOD = params['ZOD']
    if ZOD is None:
        ZOD = globals()['ZOD']

    globals()['ZOD'] += 1

    if params['color'] is None:
        pl = ax.fill_between(xdata*xscale, yscale*low_lim, yscale*up_lim, alpha=params['alpha'], linewidth=0, zorder=1)
    else:
        pl = ax.fill_between(xdata*xscale, yscale*low_lim, yscale*up_lim, facecolor=params['color'],
                             alpha=params['alpha'], linewidth=0, zorder=1)

    col = cl.rgb2hex(pl.get_facecolor()[0])

    if params['alpha_lines'] != 0:
        ax.errorbar(xdata*xscale, yscale*low_lim, color = col, linewidth=params['linewidth'], zorder = ZOD,
                    alpha = params["alpha_fill_edge"])
        ax.errorbar(xdata*xscale, yscale*up_lim, color = col, linewidth=params['linewidth'], zorder = ZOD,
                    alpha = params["alpha_fill_edge"])

    ebar = None
    if center is not None:
        ebar = ax.errorbar(xdata*xscale, center*yscale, linewidth=params['linewidth'], color=col, zorder=ZOD,
                           alpha = params['alpha_lines'], **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        if ebar is not None:
            _update_handles(ax,(ebar, pl))
        else:
            _update_handles(ax, pl)

    set_params(**params)
    if ebar is not None:
        return ebar,pl
    else:
        return pl
