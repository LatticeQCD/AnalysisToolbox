#
# plotting.py
#
# H. Sandmeyer, D. Clarke
#
# Collection of convenience tools for plotting using matplotlib.
#


import itertools, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkEqualLengths, checkType
from latqcdtools.base.utilities import isHigherDimensional, toNumpy , envector, createFilePath
from latqcdtools.base.readWrite import readTable

try:
    if 'DISPLAY' not in os.environ:
        # Check if there's no display (headless environment)
        matplotlib.use('Agg')
    else:
        # This is to avoid possible problems with colors and themes. As far as I can tell,
        # everything works with PyQt5. 
        matplotlib.use('Qt5Agg') 
except ModuleNotFoundError:
    pass


ZOD        = 10     # Global, Orders different layers 
INITIALIZE = True   # Global, flag to ensure we only initialize once
LEGEND     = False  # Global, flag to apply legend attributes when we have one
TWINX      = None   # Global, legend info is accumulated in this axis 
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

    'ax': plt,    # Axis object that is used for the plots. Default is matplotlib.pyplot.
    'ZOD': None,  # Controls where in foreground/background data/lines/bands appear.

    # Basic options affecting most plots.
    'xlabel': None,
    'alpha_xlabel': 0,           # Transparency for x-label
    'ylabel': None,
    'alpha_ylabel': 0,           # Transparency for y-label
    'xmin': None,
    'xmax': None,
    'ymin': None,
    'ymax': None,
    'title': None,
    'label': None,               # What are the data called? (Will appear in legend.)
    'color': None,               # Color for your data. (By default each new set automatically gets different color.)
    'facecolor' : 'color',       # Sometimes one differentiates between colors of faces and edges. (By default same as color.)
    'font_size': 12,             # Controls font size for everything. 
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
    'hatch': None,               # Fill pattern
    'linewidth': 1,              # Linewidth of line plots
    'capsize': 1.5,              # Length of caps af error bars
    'elinewidth': 1,             # Linewidth of the error bars of caps af error bars
    'grid' : False,              # Do you want to put a grid in the background?

    # Data markers
    'marker': "iter",            # Symbol used for plotting data. (Set to 'None' if you don't want any.)
    'markersize': 8,             # Size of the symbols.
    'markerfill': False,         # If False, markers are hollow

    # Options for the legend. 
    # legendpos:
    #   2     9     1           'upper left'    'upper center'  'upper right'
    #   6     10    7           'center left'   'center'        'center right'
    #   3     8     4           'lower left'    'lower center'  'lower right'
    'legendpos': 'best',
    'bbox_to_anchor': None,      # Manual position of the legend. If this is not None, then legendpos
                                 #   controls the point on the legend that bbox_to_anchor controls.
                                 #   for example legendpos='lower left' along with bbox_to_anchor=(0,0)
                                 #   puts the lower left corner of the legend at the bottom left corner
                                 #   of the plot. The same legendpos with bbox_to_anchor=(1,1) puts the
                                 #   lower left corner of the legend in the upper right corner of the plot.
                                 #   This is not my fault--this was a matplotlib design choice.
    'legend_ncol': 1,            # Number of columns in the legend.
    'legend_col_spacing': None,  # Spacing between columns in the legend.
    'legendframeon' : True,      # Thin border around the legend.
    'legend_title': None,        # Title of the legend.
    'handletextpad': 0.2,        # Spacing between symbol and text in legend.
    'alpha_legend': 1,           # Transperancy for the legend.

    # Adjust special aspects of the plot's axes
    'xtick_freq': None,
    'ytick_freq': None,
    'xtick_every_n': 2,    # If xtick_freq or ytick_freq is not None, label every nth tick.
    'ytick_every_n': 2,
    'xscale': 1.0,         # Scale data in xdata by this factor.
    'yscale': 1.0,
    'xlogscale': False,    # Should we use a log scale for the x-axis?
    'ylogscale': False,
    'xlogbase': 10,
    'ylogbase': 10,
    'orientation' : 'vertical',
}


# Used for later checks that the user did not pass a wrong parameter by mistake.
allowed_params = set(default_params.keys())
allowed_params = allowed_params | {'linestyle','ha','va'}


# ---------------------------------------------------------------------------------------------- SOME EXTERNAL FUNCTIONS


def latexify(bold=False):
    """ 
    Allows use of LaTeX symbols in plots. The physics package is included, allowing use of
    convenient functions like ev. 
    """
    try:
        logger.info("Using LaTeX to make pretty plots.")
        if bold:
            plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{physics}\boldmath"
        else:
            plt.rcParams['text.latex.preamble'] = r"\usepackage{lmodern}\usepackage{amssymb}\usepackage{physics}"
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'cmr10'
        plt.rcParams['axes.formatter.use_mathtext'] = True
    except Exception as e:
        logger.warn("You seem to be missing some LaTeX packages, so latexify won't work. Encountered exception:")
        logger.warn(e)


def resetLEGEND():
    global LEGEND, TWINX
    LEGEND = False
    TWINX = None


def clearPlot():
    """ 
    Clears plot object, legend handles, and zorder. Useful if you want to do multiple plots in the same script. 
    """
    logger.debug('Reset plot defaults.')
    global INITIALIZE, ZOD, LEGEND, legend_labels, legend_handles
    ZOD = 10
    INITIALIZE = True
    resetLEGEND() 
    legend_handles = { plt : [] }
    legend_labels  = { plt : [] }
    plt.clf()
    plt.close('all')


def getColorGradient(NUM_COLORS,map='viridis') -> list:
    """ 
    Generate perceptually uniform set of colors. Useful when you need more than 8 colors.

    Args:
        NUM_COLORS (int): number of colors you need 
        map (str, optional): use custom colormap. Defaults to 'viridis'.

    Returns:
        list: colors 
    """
    checkType('int',NUM_COLORS=NUM_COLORS)
    checkType(str,map=map)
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
    """ 
    Lets the user adjust the default parameter settings. 
    """
    for key, val in kwargs.items():
        default_params[key] = val


def getSubplots(x,y):
    """ 
    Get fig and axs objects when you want x figures across and y figures vertically.
    I wrapped this because matplotlib's convention for x and y is the opposite as their
    convention for figsize, which is so incredibly confusing. 
    
    Args:
        x (int): number of panels in x-direction 
        y (int): number of panels in y-direction 

    Returns:
        fig, axs: fig object, list (if 1-d) of ax objects or tuple (if 2-d)
    """
    checkType('int',x=x)
    checkType('int',y=y)
    fig, axs = plt.subplots(y,x,figsize=(4*x,4*y))
    return fig, axs


def getTwinAxes():
    """
    Get two axes objects ax1 and ax2 that share the same x-axis. Useful for graphing two functions
    on the same plot with different scales.
    """
    global TWINX
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    TWINX = ax2
    return ax1, ax2


# ---------------------------------------------------------------------------------------------- SOME INTERNAL FUNCTIONS


def _initializePlt(params):
    """ 
    Set up inital plot parameters, like its size. I tried to introduce a global variable INITIALIZE that checks
    so that this only gets called once per plot. 
    """
    checkType(dict,params=params)
    fill_param_dict(params)
    global INITIALIZE
    if INITIALIZE:
        logger.debug("Plot initializer called!")
        INITIALIZE = False
        plt.rcParams['savefig.bbox']      = 'standard'
        plt.rcParams['axes.titlesize']    = params['font_size']
        plt.rcParams['font.size']         = params['font_size']
        plt.rcParams['ytick.labelsize']   = params['font_size']
        plt.rcParams['xtick.labelsize']   = params['font_size']
        plt.rcParams['axes.labelsize']    = params['font_size']
        plt.rcParams['legend.fontsize']   = params['font_size']
        plt.rcParams['font.weight']       = params['font_weight']
        set_xrange(params['xmin'],params['xmax'])
        set_yrange(params['ymin'],params['ymax'])


def _update_labels(ax,label):
    checkType(str,label=label)
    global legend_labels, TWINX
    if ax in legend_labels:
        legend_labels[ax].append(label)
    else:
        legend_labels[ax] = [label]
    if (TWINX is not None) and (ax != TWINX): 
        if TWINX in legend_labels:
            legend_labels[TWINX].append(label)
        else:
            legend_labels[TWINX] = [label]


def _update_handles(ax,handle):
    global legend_handles, TWINX
    if ax in legend_handles:
        legend_handles[ax].append(handle)
    else:
        legend_handles[ax] = [handle]
    if (TWINX is not None) and (ax != TWINX): 
        if TWINX in legend_handles:
            legend_handles[TWINX].append(handle)
        else:
            legend_handles[TWINX] = [handle]


def _prepare_legend(ax,params):
    """
    Calls set_params, which adds labels and prepares the plot for display.
    """
    global TWINX
    if TWINX is None:
        set_params(**params)
    else:
        if ax==TWINX:
            set_params(**params)


def _getAxObject(params):
    if params['ax']==plt:
        return plt.gca()
    else:
        return params['ax']


def _getZOD(params):
    ZOD = params['ZOD']
    if ZOD is None:
        ZOD = globals()['ZOD']
    return ZOD


def _set_xmin(ax,x_min=None):
    if x_min is not None:
        try:
            x1, x2 = ax.get_xlim()
            ax.set_xlim([x_min,x2])
        except AttributeError:
            # The issue here is that the get_xlim() won't know what x2
            # is without having seen the data first.
            logger.TBRaise('Must set x/y min/max after plotting data.',exception=AttributeError)


def _set_xmax(ax,x_max=None):
    if x_max is not None:
        try:
            x1, x2 = ax.get_xlim()
            ax.set_xlim([x1,x_max])
        except AttributeError:
            logger.TBRaise('Must set x/y min/max after plotting data.',exception=AttributeError)


def _set_ymin(ax,y_min=None):
    if y_min is not None:
        try:
            y1, y2 = ax.get_ylim()
            ax.set_ylim([y_min,y2])
        except AttributeError:
            logger.TBRaise('Must set x/y min/max after plotting data.',exception=AttributeError)


def _set_ymax(ax,y_max=None):
    if y_max is not None:
        try:
            y1, y2 = ax.get_ylim()
            ax.set_ylim([y1,y_max])
        except AttributeError:
            logger.TBRaise('Must set x/y min/max after plotting data.',exception=AttributeError)


def _rescale(scale,data):
    """ 
    Rescale the data, taking into account there may be no data.

    Args:
        scale (float)
        data (array-like or None)
    """
    if data is not None:
        return np.copy(data*scale)


def fill_param_dict(params):
    """ 
    Collection of default parameters for plotting routines. If a key does not exist in params, it is defined with
    a default value.
    
    Args:
        params (dict): All parameters that are already set by the user
    """
    global LEGEND

    if not LEGEND:
        # When filling params, check if we show the legend. This is triggered by one of these keys
        # being different from the default value. 
        for key in ('legend_title', 'legendpos', 'legend_ncol', 'legendpos_col_spacing', 'label', 
                    'alpha_legend', 'bbox_to_anchor', 'legendframeon'):
            if key in params:
                if params[key] != default_params[key]:
                    logger.debug('Found legend trigger',key)
                    LEGEND = True

    for key, val in default_params.items():
        params.setdefault(key,val)


def _add_optional(params) -> dict:
    """ 
    Optional parameters not included in _fill_param_dict.

    Args:
        **params: Additional parameters that can be set.

    Returns:
        dict: also optional parameters 
    """
    reference = {}
    fill_param_dict(reference)
    ret = {}
    for key in params:
        if not key in reference:
            ret[key] = params[key]
    return ret


def set_params(**params):
    """ 
    Set additional parameters to the plot. For example set a title or label.

    Args:
        **params: Additional parameters that can be set.
    """
    global LEGEND

    fill_param_dict(params)
    ax = _getAxObject(params)

    # Check for some contradictory settings.
    if (params['xlabelpos'] is not None) and (params['xlabel'] is None):
        logger.warn('Set xlabelpos with no xlabel.')
    if (params['ylabelpos'] is not None) and (params['ylabel'] is None):
        logger.warn('Set ylabelpos with no ylabel.')
    if params['xlogscale'] and params['xtick_freq'] is not None:
        logger.warn("xtick_freq assumes no log scale.")
    if params['ylogscale'] and params['ytick_freq'] is not None:
        logger.warn("ytick_freq assumes no log scale.")

    if params['xlabel'] is not None:
        checkType(str,xlabel=params['xlabel'])
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
        checkType(str,ylabel=params['ylabel'])
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

    ax.tick_params(axis='both', which='major', labelsize=params['font_size']) 

    if params['title'] is not None:
        checkType(str,title=params['title'])
        plt.title(params['title'])

    if params['xlogscale']:
        ax.set_xscale('log',base=params['xlogbase'])

    if params['ylogscale']:
        ax.set_yscale('log',base=params['ylogbase'])

    if params['ticksintoplot']:
        ax.tick_params(which='both',direction='in')

    if params['surroundWithTicks']:
        ax.tick_params(which='minor',left=True,bottom=True,top=True,right=True)
        ax.tick_params(which='major',left=True,bottom=True,top=True,right=True)

    if params['grid']:
        ax.grid(True, color='lightgray',linestyle='dotted',zorder=BACKGROUND)

    set_xrange(params['xmin'],params['xmax'],ax)
    set_yrange(params['ymin'],params['ymax'],ax)

    if LEGEND:
        if not ax in legend_handles:
            logger.TBRaise('Legend for axis',ax,'was activated without any label.')
        leg = ax.legend(legend_handles[ax], legend_labels[ax], numpoints=1, bbox_to_anchor = params['bbox_to_anchor'],
                        title=params['legend_title'], loc=params['legendpos'], ncol=params['legend_ncol'],
                        columnspacing=params['legend_col_spacing'],handletextpad = params['handletextpad'],
                        frameon = params['legendframeon'])
        leg.get_frame().set_alpha(params['alpha_legend'])
        leg.set_zorder(FOREGROUND)

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


def preliminary(x,y,text='PRELIMINARY',**kwargs):
    """ 
    Generate a PRELIMINARY tag on the plot.

    Args:
        x (float): x-position of bottom-left corner (in units of x-axis) 
        y (float): y-position of bottom-left corner (in units of y-axis) 
        text (str, optional): Text indicating result is preliminary. Defaults to 'PRELIMINARY'.
    """
    checkType(str,text=text)
    _initializePlt(kwargs)
    if 'color' in kwargs: 
        color=kwargs['color']
    else:
        color='gray'
    axObj=plt
    if 'ax' in kwargs:
        axObj=kwargs['ax']
    axObj.text(x, y, text, color=color)


def plot_file(filename, xcol=0, ycol=1, yecol=None, xecol=None, func = None, func_args = (), style='dots', **params):
    """ 
    Plot data in file. You can set the style with the style argument. Columns indexed from 0.

    Args:
        filename (str)
        xcol (int, optional): Which column is xdata. Defaults to 1.
        ycol (int, optional): Which column is ydata. Defaults to 2.
        yecol (int, optional): Which column has y error. Defaults to None.
        xecol (int, optional): Which column has x error. Defaults to None.
        func (function, optional): Apply this function to the data. Defaults to None.
        func_args (tuple, optional): Arguments to func. Defaults to ().
        style (str, optional): Choose from dots, lines, fill, and band. Defaults to 'dots'.
        **params: Additional parameters that can be set.
    """
    checkType('int',xcol=xcol)
    checkType('int',ycol=ycol)
    _initializePlt(params)
    data   = readTable(filename,dtype=str)
    xdata  = data[xcol].astype(float)
    ydata  = data[ycol].astype(float)
    yedata = None
    xedata = None

    if yecol is not None:
        yedata = data[yecol].astype(float)
    else:
        if style == "fill":
            logger.TBRaise("Need error column for filled plotting")
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
        plot_dots(xdata, ydata, yedata=yedata, xedata=xedata, **params)
    elif style == "lines":
        plot_lines(xdata, ydata, yedata=yedata, xedata=xedata, **params)
    elif style == "fill":
        plot_fill(xdata, ydata, yedata=yedata, **params)
    else:
        logger.TBRaise("Unknown style",style)


def plot_vspan(minVal,maxVal,**params):
    """ 
    Plot a vertical band.

    Args:
        minVal (float)
        maxVal (float)
    """
    _initializePlt(params)
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)
    handle   = ax.axvspan(xmin=minVal, xmax=maxVal,color=params['color'],alpha=params['alpha'],zorder=ZOD,**optional)
    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,handle)
        _prepare_legend(ax,params)
    globals()['ZOD'] += 1


def plot_hspan(minVal,maxVal,**params):
    """ 
    Plot a horizontal band.

    Args:
        minVal (float)
        maxVal (float)
    """
    _initializePlt(params)
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)
    handle   = ax.axhspan(ymin=minVal, ymax=maxVal,color=params['color'],alpha=params['alpha'],zorder=ZOD,**optional)
    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,handle)
        _prepare_legend(ax,params)
    globals()['ZOD'] += 1


def plot_hline(y,minVal=None,maxVal=None,**params):
    """ 
    Plot a horizontal line at y. 
    """
    _initializePlt(params)
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)
    x1, x2   = ax.get_xlim()
    if minVal is not None:
        x1=minVal
    if maxVal is not None:
        x2=maxVal
    logger.debug('x1, x2 = ',x1,x2)
    # The matplotlib people are psychopaths. I couldn't get this to work in any other
    # way. Yes, the comma after handle in the else branch is mandatory. 
    if (minVal is None) and (maxVal is None):
        handle = ax.axhline(y=y,color=params['color'],zorder=ZOD,alpha=params['alpha'], **optional) 
    else:
        handle, = ax.plot([x1, x2], [y, y], color=params['color'],zorder=ZOD,alpha=params['alpha'], **optional) 
    globals()['ZOD'] += 1
    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,handle)
        _prepare_legend(ax,params)


def plot_vline(x,minVal=None,maxVal=None,**params):
    """ 
    Plot a vertical line at x. 
    """
    _initializePlt(params)
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)
    y1, y2   = ax.get_ylim()
    if minVal is not None:
        y1=minVal
    if maxVal is not None:
        y2=maxVal
    logger.debug('y1, y2 = ',y1,y2)
    if (minVal is None) and (maxVal is None):
        handle = ax.axvline(x=x,color=params['color'],zorder=ZOD,alpha=params['alpha'], **optional) 
    else:
        handle, = ax.plot([x, x], [y1, y2], color=params['color'], zorder=ZOD, alpha=params['alpha'], **optional) 
    globals()['ZOD'] += 1
    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,handle)
        _prepare_legend(ax,params)


def plot_dots(xdata, ydata, yedata = None, xedata = None, **params):
    """ 
    Plot ydata vs xdata as dots. 

    Args:
        xdata (array-like)
        ydata (array-like)
        yedata (array-like, optional): y error. Defaults to None.
        xedata (array-like, optional): x error. Defaults to None.
        **params: Additional parameters that can be set.
    """

    if len(envector(yedata)) == 2:
        checkEqualLengths(xdata,ydata,yedata.T)
    else:
        checkEqualLengths(xdata,ydata,yedata)

    if len(envector(xedata)) == 2:
        checkEqualLengths(xdata,ydata,xedata.T)
    else:
        checkEqualLengths(xdata,ydata,xedata)

    xdata, ydata, xedata, yedata = toNumpy(xdata, ydata, xedata, yedata)
    _initializePlt(params)
    optional = _add_optional(params)

    # If you want lines between your data you should be using plot_lines.
    if 'linestyle' in optional:
        logger.warn('Ignoring linestyle',optional['linestyle'])
        del optional['linestyle']

    ax     = _getAxObject(params)
    xedata = _rescale(abs(params['xscale']),xedata)
    yedata = _rescale(abs(params['yscale']),yedata)
    ZOD    = _getZOD(params)

    marker = params['marker']
    if marker == "iter":
        marker = next(markers)
    if params['markerfill']:
        markerfill=params['color']
    else:
        markerfill="None"

    if params['color'] is not None:
        handle = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], yerr=yedata, xerr=xedata, marker=marker, 
                           linestyle='None', linewidth=params['linewidth'], alpha=params['alpha'], color=params['color'],
                           zorder=ZOD, markersize=params['markersize'], capsize=params['capsize'],
                           elinewidth=params['elinewidth'],
                           markerfacecolor = markerfill, **optional)
    else:
        handle = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], yerr=yedata, xerr=xedata, marker=marker, 
                           linestyle='None', linewidth=params['linewidth'], alpha=params['alpha'], zorder=ZOD,
                           markersize=params['markersize'], capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markerfacecolor = markerfill, **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,handle)
        _prepare_legend(ax,params)

    globals()['ZOD'] += 1


def plot_bar(xdata, ydata, width=None, align='edge', edgecolor='#666677',linewidth=0.2, **params):
    """ 
    Plot ydata vs xdata as bars.

    Args:
        xdata (array-like)
        ydata (array-like)
        width (float, optional): Width of bar. Defaults to xdata[1]-xdata[0].
        align (str, optional): How to align the bars. Defaults to 'edge'.
        edgecolor (str, optional): Color of bar edges. Defaults to '#666677'.
        linewidth (float, optional): Defaults to 0.2.
        **params: Additional parameters that can be set.
    """
    checkEqualLengths(xdata,ydata)
    xdata, ydata = toNumpy(xdata, ydata)
    _initializePlt(params) 
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)

    if width is None:
        width = xdata[1] - xdata[0]

    if params['color'] is not None:
        bar = ax.bar(xdata, ydata, color=params['color'], zorder=ZOD, width=width, align=align, edgecolor=edgecolor,
                     linewidth=linewidth, alpha=params['alpha'], **optional)
    else:
        bar = ax.bar(xdata, ydata, zorder=ZOD, width=width, align=align, edgecolor=edgecolor,
                     linewidth=linewidth, alpha=params['alpha'], **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,bar)
        _prepare_legend(ax,params)
    globals()['ZOD'] += 1


def plot_hist(data, bins = None, density=False, label=None, weights=None, **params):
    """ 
    Create a histogram of the array data. If you would like to plot multiple data sets in the same histogram,
    simply pass as a list or tuple of arrays of data, like data = [list1, list2, ...].

    Args:
        data (array-like)
        bins (int, optional): Number of bins. Defaults to None, which sets the number of bins automatically.
        **params: Additional parameters that can be set.
    """
    _initializePlt(params) 
    ax = _getAxObject(params)
    if bins is None:
        bins = 'auto'
        if weights is not None:
            logger.TBRaise('Matplotlib does not support automatic bins with weights')
    if isHigherDimensional(data):
        ax.hist(data, bins=bins, density=density, label=label, orientation=params['orientation'], alpha=params['alpha'],
                weights=weights)
    else:
        ax.hist(data, bins=bins, density=density, label=label, color=params['color'], orientation=params['orientation'],
                alpha=params['alpha'], weights=weights)
    if density:
        ax.set_yticklabels([])
    if label is not None:
        ax.legend()


def plot_lines(xdata, ydata, yedata=None, xedata=None, **params):
    """ 
    Plot ydata vs xdata using lines.

    Args:
        xdata (array-like)
        ydata (array-like)
        yedata (array-like, optional): y error. Defaults to None.
        xedata (array-like, optional): x error. Defaults to None.
        **params: Additional parameters that can be set.
    """
    checkEqualLengths(xdata,ydata,yedata,xedata)
    xdata, ydata, xedata, yedata = toNumpy(xdata, ydata, xedata, yedata)
    _initializePlt(params)
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)

    xedata = _rescale(params['xscale'],xedata)
    yedata = _rescale(params['yscale'],yedata)

    marker = params['marker']
    if marker == "iter":
        marker = next(markers)
    if params['markerfill']:
        markerfill=params['color']
    else:
        markerfill="None"

    if params['color'] is not None:
        handle = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], yerr=yedata, xerr=xedata, marker=marker, 
                           linestyle='None', linewidth=params['linewidth'], color=params['color'], zorder=ZOD,
                           markersize=params['markersize'], capsize=params['capsize'], elinewidth=params['elinewidth'],
                           markeredgewidth=params['elinewidth'], alpha=params['alpha_dots'], 
                           markerfacecolor = markerfill)
    else:
        handle = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], yerr=yedata, xerr=xedata, marker=marker, 
                           linestyle='None', linewidth=params['linewidth'], zorder=ZOD, markersize=params['markersize'],
                           capsize=params['capsize'], elinewidth=params['elinewidth'], 
                           markeredgewidth=params['elinewidth'], alpha=params['alpha_dots'],
                           markerfacecolor = markerfill)

    globals()['ZOD'] += 1

    col = handle[0].get_color()

    line = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], color = col, linewidth=params['linewidth'], 
                       zorder = ZOD, alpha = params["alpha"], **optional)

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,(line,handle))
        _prepare_legend(ax,params)


def plot_fill(xdata, ydata, yedata, xedata=None, center=True, **params):
    """ 
    Plot a filled region within ydata +/- yedata. Can set xedata along with yedata=None for vertical bands.

    Args:
        xdata (array-like)
        ydata (array-like)
        yedata (array-like): y error 
        xedata (array-like, optional): x error. Defaults to None.
        center (bool): Do you show the central line? Defaults to True. 
        **params: Additional parameters that can be set.
    """
    if (yedata is None) and (xedata is None):
        logger.TBRaise("Please pass some error bars.")
    if (yedata is not None) and (xedata is not None):
        logger.TBRaise("Please pass either x-error or y-error, not both.")
    xdata, ydata, xedata, yedata = toNumpy(xdata, ydata, xedata, yedata)

    if len(envector(yedata)) == 2:
        checkEqualLengths(xdata,ydata,xedata,yedata.T)
    else:
        checkEqualLengths(xdata,ydata,xedata,yedata)

    _initializePlt(params) 
    optional = _add_optional(params)
    ax       = _getAxObject(params)
    ZOD      = _getZOD(params)

    if center :
        if params['color'] is not None:
            handle = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], linewidth=params['linewidth'], 
                                 color=params['color'], zorder=ZOD+1, alpha = params['alpha_lines'], **optional)
        else:
            handle = ax.errorbar(xdata*params['xscale'], ydata*params['yscale'], linewidth=params['linewidth'], 
                                 zorder=ZOD+1, alpha = params['alpha_lines'], **optional)
        col = handle[0].get_color()
    else:
        col = params['color']

    if params['facecolor'] == 'color':
        facecol=col
    else:
        facecol=params['facecolor']

    if xedata is None:
        if len(yedata) == 2:
            pl = ax.fill_between(xdata*params['xscale'],
                                 (np.asarray(ydata*params['yscale']) - np.asarray(yedata[0]*params['yscale'])),
                                 (np.asarray(ydata*params['yscale']) + np.asarray(yedata[1]*params['yscale'])), 
                                 facecolor=facecol, alpha=params['alpha'], linewidth=params['linewidth'], 
                                 zorder=ZOD, edgecolor=col, hatch=params['hatch'])
        else:    
            pl = ax.fill_between(xdata*params['xscale'],
                                 (np.asarray(ydata*params['yscale']) - np.asarray(yedata*params['yscale'])),
                                 (np.asarray(ydata*params['yscale']) + np.asarray(yedata*params['yscale'])), 
                                 facecolor=facecol, alpha=params['alpha'], linewidth=params['linewidth'], 
                                 zorder=ZOD, edgecolor=col, hatch=params['hatch'])
    else:
        pl = ax.fill_betweenx(ydata*params['yscale'], 
                              (np.asarray(xdata*params['xscale']) - np.asarray(xedata*params['xscale'])),
                              (np.asarray(xdata*params['xscale']) + np.asarray(xedata*params['xscale'])), 
                              facecolor=facecol, alpha=params['alpha'],linewidth=params['linewidth'], 
                              zorder=ZOD, edgecolor=col, hatch=params['hatch'])
    globals()['ZOD'] += 2

    if params['label'] is not None:
        _update_labels(ax,params['label'])
        _update_handles(ax,(handle,pl))
        _prepare_legend(ax,params)


def plot_matrix(mat,vmin=None,vmax=None):
    """ Plot matrix as a heatmap.

    Args:
        mat (np.ndarray): correlation matrix
        ax (matplotlib ax object): Defaults to plt.
    """
    checkType(np.ndarray,mat=mat)
    heatmapColors = ['blue', 'white', 'red']
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom_cmap', heatmapColors)
    plt.imshow(mat, cmap=cmap, vmin=vmin,vmax=vmax)
    plt.colorbar()


def saveFigure(filename,**kwargs):
    """
    Wrapper for plt.savefig that creates the directory path if it doesn't exist already.
    """
    createFilePath(filename)
    plt.savefig(filename,**kwargs)
