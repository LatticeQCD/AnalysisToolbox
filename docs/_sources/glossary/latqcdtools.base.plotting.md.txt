latqcdtools.base.plotting
=============

```Python
_add_optional(params) -> dict:
'''
Optional parameters not included in _fill_param_dict.

Args:
    **params: Additional parameters that can be set.

Returns:
    dict: also optional parameters 
'''
```
```Python
_getAxObject(params):
'''
'''
```
```Python
_initializePlt(params):
'''
Set up inital plot parameters, like its size. I tried to introduce a global variable INITIALIZE that checks
so that this only gets called once per plot. 
'''
```
```Python
_rescale(scale, data):
'''
Rescale the data, taking into account there may be no data.

Args:
    scale (float)
    data (array-like or None)
'''
```
```Python
_set_xmax(ax, x_max=None):
'''
'''
```
```Python
_set_xmin(ax, x_min=None):
'''
'''
```
```Python
_set_ymax(ax, y_max=None):
'''
'''
```
```Python
_set_ymin(ax, y_min=None):
'''
'''
```
```Python
_update_handles(ax, handle):
'''
'''
```
```Python
_update_labels(ax, label):
'''
'''
```
```Python
clearPlot():
'''
Clears plot object, legend handles, and zorder. Useful if you want to do multiple plots in the same script. 
'''
```
```Python
fill_param_dict(params):
'''
Collection of default parameters for plotting routines. If a key does not exist in params, it is defined with
a default value.

Args:
    params (dict): All parameters that are already set by the user
'''
```
```Python
getColorGradient(NUM_COLORS, map='viridis') -> list:
'''
Generate perceptually uniform set of colors. Useful when you need more than 8 colors.

Args:
    NUM_COLORS (int): number of colors you need 
    map (str, optional): use custom colormap. Defaults to 'viridis'.

Returns:
    list: colors 
'''
```
```Python
getSubplots(x, y):
'''
Get fig and axs objects when you want x figures across and y figures vertically.
I wrapped this because matplotlib's convention for x and y is the opposite as their
convention for figsize, which is so incredibly confusing. 

Args:
    x (int): number of panels in x-direction 
    y (int): number of panels in y-direction 

Returns:
    fig, axs: fig object, list (if 1-d) of ax objects or tuple (if 2-d)
'''
```
```Python
latexify(bold=False):
'''
Allows use of LaTeX symbols in plots. The physics package is included, allowing use of
convenient functions like ev. 
'''
```
```Python
plot_bar(xdata, ydata, width=None, align='edge', edgecolor='#666677', linewidth=0.2, **params):
'''
Plot ydata vs xdata as bars.

Args:
    xdata (array-like)
    ydata (array-like)
    width (float, optional): Width of bar. Defaults to xdata[1]-xdata[0].
    align (str, optional): How to align the bars. Defaults to 'edge'.
    edgecolor (str, optional): Color of bar edges. Defaults to '#666677'.
    linewidth (float, optional): Defaults to 0.2.
    **params: Additional parameters that can be set.
'''
```
```Python
plot_dots(xdata, ydata, yedata=None, xedata=None, **params):
'''
Plot ydata vs xdata as dots. 

Args:
    xdata (array-like)
    ydata (array-like)
    yedata (array-like, optional): y error. Defaults to None.
    xedata (array-like, optional): x error. Defaults to None.
    **params: Additional parameters that can be set.
'''
```
```Python
plot_file(filename, xcol=0, ycol=1, yecol=None, xecol=None, func=None, func_args=(), style='dots', **params):
'''
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
'''
```
```Python
plot_fill(xdata, ydata, yedata, xedata=None, center=True, **params):
'''
Plot a filled region within ydata +/- yedata. Can set xedata along with yedata=None for vertical bands.

Args:
    xdata (array-like)
    ydata (array-like)
    yedata (array-like): y error 
    xedata (array-like, optional): x error. Defaults to None.
    center (bool): Do you show the central line? Defaults to True. 
    **params: Additional parameters that can be set.
'''
```
```Python
plot_hist(data, bins=None, density=False, label=None, weights=None, **params):
'''
Create a histogram of the array data. If you would like to plot multiple data sets in the same histogram,
simply pass as a list or tuple of arrays of data, like data = [list1, list2, ...].

Args:
    data (array-like)
    bins (int, optional): Number of bins. Defaults to None, which sets the number of bins automatically.
    **params: Additional parameters that can be set.
'''
```
```Python
plot_hline(y, minVal=None, maxVal=None, **params):
'''
Plot a horizontal line at y. 
'''
```
```Python
plot_hspan(minVal, maxVal, **params):
'''
Plot a horizontal band.

Args:
    minVal (float)
    maxVal (float)
'''
```
```Python
plot_lines(xdata, ydata, yedata=None, xedata=None, **params):
'''
Plot ydata vs xdata using lines.

Args:
    xdata (array-like)
    ydata (array-like)
    yedata (array-like, optional): y error. Defaults to None.
    xedata (array-like, optional): x error. Defaults to None.
    **params: Additional parameters that can be set.
'''
```
```Python
plot_matrix(mat, vmin=None, vmax=None):
'''Plot matrix as a heatmap.

Args:
    mat (np.ndarray): correlation matrix
    ax (matplotlib ax object): Defaults to plt.
'''
```
```Python
plot_vline(x, minVal=None, maxVal=None, **params):
'''
Plot a vertical line at x. 
'''
```
```Python
plot_vspan(minVal, maxVal, **params):
'''
Plot a vertical band.

Args:
    minVal (float)
    maxVal (float)
'''
```
```Python
preliminary(x, y, text='PRELIMINARY', **kwargs):
'''
Generate a PRELIMINARY tag on the plot.

Args:
    x (float): x-position of bottom-left corner (in units of x-axis) 
    y (float): y-position of bottom-left corner (in units of y-axis) 
    text (str, optional): Text indicating result is preliminary. Defaults to 'PRELIMINARY'.
'''
```
```Python
resetLEGEND():
'''
'''
```
```Python
saveFigure(filename, **kwargs):
'''
Wrapper for plt.savefig that creates the directory path if it doesn't exist already.
'''
```
```Python
set_default_param(**kwargs):
'''
Lets the user adjust the default parameter settings. 
'''
```
```Python
set_params(**params):
'''
Set additional parameters to the plot. For example set a title or label.

Args:
    **params: Additional parameters that can be set.
'''
```
```Python
set_xrange(xmin=None, xmax=None, ax=<module 'matplotlib.pyplot' from '/home/dclarke/.local/lib/python3.13/site-packages/matplotlib/pyplot.py'>):
'''
'''
```
```Python
set_yrange(ymin=None, ymax=None, ax=<module 'matplotlib.pyplot' from '/home/dclarke/.local/lib/python3.13/site-packages/matplotlib/pyplot.py'>):
'''
'''
```
