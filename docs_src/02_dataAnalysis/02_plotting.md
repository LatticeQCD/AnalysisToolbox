# Plotting 

## Making Plots

The `analysistoolbox` also has some wrappers for Matplotlib, which allow you to make Matplotlib 
plots easily using python. All of them can be found in 
```Python
import latqcdtools.base.plotting
```
It is recommended to use this module because you can add all the power and flexibility of python 
to your plotting routines.

## The basics

You probably want to get started quickly, so here's a minimal example. **Please be mindful of the 
order of the commands.**
```Python
latexify()
plot_file('susceptibility.d',1,2,3,label='$56^3\\times8~~m_s/m_l=80$',color=colors[1],marker=markers_1[1])
set_params(xlabel="T [MeV]",
           ylabel="$N_s^3\\,\\chi^{\\rm{bare}}$",
           title ="$N_\\tau=8$")
plt.savefig("susceptibility.pdf")
plt.show()
```
The above python code will create a plot using the data in `susceptibility.d`, save it in 
`susceptibility.pdf`, and also show the plot to you. Let us go through each line. You probably 
noticed that the labels and titles use $\LaTeX$ syntax. In order for this to work, a call to
```Python
latexify()
```
is needed at the beginning. In addition to using nice $\LaTeX$ font this command also automatically 
puts the axes in scientific notation if that would look better. I recommend that you always use this 
command in your plots. The command
```Python
set_params()
```
lets you specify things needed for every plot such as titles and axis labels. The `set_params` 
command should always come _after_ your `plot_file` calls and _before_ your `plt.savefig` and 
`plt.show` calls. Somehow the `plot_file` seems to overwrite the `set_params`. For a full list 
of parameters have a look at this module. The command
```Python
plot_file(filename, xcol, ycol, xecol, yecol=None)
```
is what adds the data to the plot. The arguments `xcol`, `ycol`, `xecol`, and `yecol` are the column 
numbers of the data file for the x-data, y-data, x error data, and y error data, respectively. After 
this there follows some optional arguments. The `marker` is the symbol used to mark the data on the 
plot, and the `color` is the color of those markers. There are 8 colors in the `colors` list; 
other possible colors to use are given in the module. There are 11 markers in the `markers` list, 
and since this module is based on Matplotlib, you find other possibilities for markers 
[here](https://matplotlib.org/api/markers_api.html). The last two commands are self-explanatory. 
In case you wish to change the style of your plot, this can be accomplished with the
```Python
style='dots'
```
option in your `plot_file` arguments. The allowed styles include `dots`, `lines`, `fill`, and 
`band`. **Just as a friendly tip: If you use a style other than `dots`, which is the default, 
you may need to use `marker=None` or `marker=','` to suppress the data points. This will be 
necessary if the data points are close enough together and larger than the line itself.**

Again since we want to be able to plot quickly within python scripts, you might find yourself in 
a situation where you want to plot list elements rather than a file. As a minimal example, one can 
use
```Python
plot_dots(xdata,ydata,yedata=edata)
```
instead of `plot_file` to plot the list `ydata` with its error bars `edata` against `xdata`. 
Other plot style options besides `plot_dots` include `plot_lines`, `plot_fill`, and `plot_band`.

If you are using a single python script to plot multiple files, python will remember the labels 
from the last plot for some reason. Thankfully there is a command which addresses this problem: 
After you have written your code to generate some plot `plot1`, you can call `clear_legend_labels()` 
before your code to generate `plot2`. This will reset the legend and generate your correct 
plot as expected.

## Some options you might commonly require

### Changing axis ranges

The ranges of plots are chosen automatically, but this automatic choice may not always suit your 
purposes. In this can you can use the commands
```Python
set_xmin(-1)
set_xmax(1)
set_ymin(-2)
set_ymax(2)
```
to set the x-range to lie in the interval [-1,1] and the y-range to lie in the interval [-2,2]. 
Please note that the axis ranges are NOT controlled by `set_params`. `set_params` does have, for 
example, the parameter `xmin`, but the parameter `xmin` contained in `set_params` controls only 
what data is used for the plot, and does not control the plot range directly. Please also put 
these commands after `plot_file` but before `plt.savefig` and `plt.show`.

### Changing the location of the legend

The legend position is normally chosen automatically, but you may have other ideas. The 
`legendpos` parameter in `set_params` can be used to adjust the position of the legend. 
For example
```Python
legendpos =  5
```
will put the legend on the right, within the boundary of the plot. A full list of positions that 
lie within the boundary of the plot can be found in the plotting module. You may need to put the 
legend outside of the plot. In that case use
```Python
legendpos =  (1,0)
```
to put the legend on the right, outside the plot, at the bottom. You can increase the second 
coordinate to put the legend higher up.

### Vertical and horizontal lines; bands

Sometimes you may wish to add a vertical or horizontal line to indicate the position of some 
important quantity. These can be created with, for example
```Python
plt.axvline(x=0, color=colors[0])
plt.axhline(y=0, color=colors[0])
```
You also have the option of making them dashed using `linestyle="dashed"`. You may also want 
to create an error band or something like this. Matplotlib also has a nice command for this, 
namely
```Python
plt.axvspan(xmin, xmax, alpha=0.5, color=colors[0])
plt.axhspan(ymin, ymax, alpha=0.5, color=colors[0])
```
Here the `min` and `max` arguments are used to specify the width of the band, and 0 < `alpha` < 1 
is the opacity. If you have a file with data and error bars, and you would like to plot these as 
error bands, you can use
```Python
style='fill'
```
as an argument to your `plot_file` command.


### Changing the size of the plot

Maybe you have a computer with high resolution, the image might come up rather small when you 
`plt.show()`. To fix this, you can use `latexify` to specify the width, height, and text size of 
your figure as, for example `latexify(fig_width=40,fig_height=28,font_size=18)`.

### Printing only an inset to pdf

When you use `plt.show()` to look at your figure, a small region in the graph may catch your 
interest, and you might zoom in on that. It is possible to print this inset to pdf. 
Use
```Python
plt.show()
plt.ion()
plt.savefig("example.pdf")
```
Note that the relative order of `plt.show()` and `plt.savefig()` has changed.

### Adding a text box

Do you have a preliminary result, which requires you to write PRELIMINARY in? You can accomplish 
this with
```Python
plt.text(xcoordinate, ycoordinate, 'PRELIMINARY', color='gray')
```
The `xcoordinate` and `ycoordinate` are measured according to the axes themselves; for example if 
my x-axis is temperature and I want to put the text at 150 [MeV], I just set `xcoordinate=150`. 
Further options for `text` can be found 
[here](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html).

### Rescaling an axis

We need to do this pretty often, for example when changing from lattice units to physical units. To 
rescale the x-axis by a factor 2, pass to your `plot_` command `xscale=2`. `yscale` rescales the y-axis.
