# 
# testHist.py                                                               
# 
# D. Clarke
# 
# Quick test of histogram method. 
#

from latqcdtools.base.plotting import plt, latexify, plot_hist, colors 
from latqcdtools.interfaces.interfaces import loadGPL

latexify()

bareCorr = loadGPL('hist2.gpx')
factor   = 9*8.3724494170096e-08
data2    = bareCorr[:,6]*factor
bareCorr = loadGPL('hist1.gpx')
data1    = bareCorr[:,6]

data   = (data1,data2)
labels = ['hist1','hist2']

plot_hist(data,bins=20,density=True,xlabel='corr$_{\\rm disc}(0)$',color=colors[0],label=labels) 

plt.show()