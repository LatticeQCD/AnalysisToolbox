# 
# test_plot.py                                                               
# 
# H. Sandmeyer
# 
# Test several of the methods in the plotting.py module. This script is a good
# example script to see what kind of fanciful things you can do with plots.
#

import numpy as np
import matplotlib.pyplot as plt
from latqcdtools.base.plotting import latexify, plot_file, plot_bar, plot_func

latexify()

plot_file("data1.txt", 1, 2, 3, 4, label="gauss noise", xlabel="$\\frac{1}{z}$", capsize = 2)
plot_file("wurf.dat", 1, 3, 4, 2, style="lines", label="quadratic", marker = 'o')
plot_file("data2.txt", 1, 2, 3, 4, style="fill", alpha=0.3, label="linear")
plot_file("data3.txt", 1, 2, 3, 4, style="fill", alpha=0.9, label="linear 2", legendpos = "best", capsize = 2)
plot_bar([-6,-5],[1,2])
plot_func(np.sin)
plt.show()
plt.close()
