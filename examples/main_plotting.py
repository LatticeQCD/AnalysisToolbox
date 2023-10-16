# 
# main_plotting.py                                                               
# 
# D. Clarke 
# 
# Some explicit examples how to use the plotting module. The plotting module inherits from
# matplotlib.pyplot, so a lot of the options and methods of matplotlib can be easily folded in. 
# 

import numpy as np
from latqcdtools.base.plotting import latexify, set_params, plot_file, clearPlot, preliminary, \
    plot_bar, plot_hist, plt
from latqcdtools.base.initialize import initialize, finalize
from latqcdtools.statistics.statistics import plot_func
from latqcdtools.interfaces.interfaces import loadGPL


initialize('example_plotting.log')

latexify()

# To start we plot a file, wurf.dat. After the file name, we specify, x, y, y-err, and x-err
# columns to use from the file. For the style you can choose between dots, line, fill, band.
# You should in general be able to pass standard matplotlib keyword arguments like 'marker'
# to many of the plotting module's methods. 
plot_file("wurf.dat", 0, 2, 3, 1, style="dots", label="wurf", marker='o')

# Sometimes a result is preliminary. You can put a marker that indicates your plot is preliminary
# using this command. The arguments set the x- and y-coordinates in units of the axes. 
preliminary(x=1,y=2)

# The set_params command should in general be the last thing you do before plotting or saving
# a figure. This is because some plotting methods may change these parameters behind the scenes.
# For a full list of parameter options, check out the default_params dictionary in the
# plotting module. I would also like to point out that the latexify() command allows use of
# LaTeX macros from the physics.sty package, for instance \ev and \int.
set_params(xlabel   ="$T$ [MeV]",
           ylabel   ="$\\chi^{\\rm{bare}}_{\\Re P}$",
           title    ="$N_\\tau=8$",
           legendpos=(1,0),
           legend_title="$\\ev{\\int dx}$",
           xlabelpos=(0.6,0.05),
           ylabelpos=(0.05,0.4),
           font_size=10,
           xmin=1.0,
           xmax=3.0,
           ymin=1.0
          )

plt.show()

# If you would like to use one script to generate multiple independent plots, such as done in this
# example script, you should call clearPlot() in between each new plot.
clearPlot()

# Here we show some more plotted data, this time also adjusting transparency (alpha).
plot_file("data1.txt", 0, 1, 2, 3, label="gauss noise", xlabel="$\\frac{1}{z}$", capsize = 2)
plot_file("wurf.dat" , 0, 2, 3, 1, style="lines", label="quadratic", marker = 'o')
plot_file("data2.txt", 0, 1, 2, 3, style="fill", alpha=0.3, label="linear")
plot_file("data3.txt", 0, 1, 2, 3, style="fill", alpha=0.9, label="linear 2", legendpos = "best", capsize = 2)

# At the moment plot_func is inside the statistics.statistics module to allow for error propagation.
plot_func(np.sin,domain=(-5,10))

set_params(legendpos=9)
plt.show()
clearPlot()

# Trivial bar plot.
plot_bar([0,1,2,3,4,5,6,7,8],[32500,32000,32600,35200,36000,37655,31433,31917,36188])
set_params(ylabel='Euros')
plt.show()
clearPlot()

# Simple histogram. Note that one can read in gpl-type files like used in Peter LePage's code.
bareCorr = loadGPL('hist2.gpx')
factor   = 9*8.3724494170096e-08
data2    = bareCorr[:,6]*factor
bareCorr = loadGPL('hist1.gpx')
data1    = bareCorr[:,6]

# Here we put two histograms on the same plot.
data   = (data1,data2)
labels = ['hist1','hist2']
plot_hist(data,bins=20,density=True,xlabel='corr$_{\\rm disc}(0)$',label=labels) 

# You always have to call savefig before show, otherwise savefig prints a blank.
plt.savefig('exampleHist.pdf')

plt.show()

finalize()