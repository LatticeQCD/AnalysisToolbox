from latqcdtools.base.plotting import *

latexify()

QM=True

if QM:
  plot_file("pressure_HRG_muB0.00_b0.40_QM",xcol=1,ycol=3,marker=None,style="lines",linewidth="0.6",label="$\\mu_B=0.0$")
  plot_file("pressure_HRG_muB1.00_b0.40_QM",xcol=1,ycol=3,marker=None,style="lines",linewidth="0.6",label="$\\mu_B=1.0$")
  set_params(ylabel="$P/T^4$",xlabel="$T$ [MeV]",title="QMHRG")
else:
  plot_file("pressure_HRG_muB0.00_b0.40_QM",xcol=1,ycol=2,marker=None,style="lines",linewidth="0.6",label="$\\mu_B=0.0$")
  plot_file("pressure_HRG_muB1.00_b0.40_QM",xcol=1,ycol=2,marker=None,style="lines",linewidth="0.6",label="$\\mu_B=1.0$")
  set_params(ylabel="$P/T^4$",xlabel="$T$ [MeV]",title="PDGHRG")


plt.show()

