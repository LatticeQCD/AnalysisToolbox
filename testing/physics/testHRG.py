# 
# testHRG.py                                                               
# 
# D. Clarke, J. Goswami 
# 
# A test for some of the HRG methods. 
#

import numpy as np
import matplotlib.pyplot as plt
import latqcdtools.base.logger as logger
from latqcdtools.base.check import print_results
from latqcdtools.base.cleanData import excludeAtCol,restrictAtCol
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.plotting import plot_lines,plot_file,set_params,latexify,colors,clear_legend_labels
from latqcdtools.base.utilities import timer
from latqcdtools.base.speedify import parallel_function_eval
from latqcdtools.math.num_deriv import diff_deriv
from latqcdtools.physics.HRG import HRG,EVHRG,HRGexact


times = timer()


EPSILON = 1e-6
SHOW_PLOTS = False # In case you want a visual to see how comparisons with older results look.

if SHOW_PLOTS:
    latexify()


T = np.linspace(130, 179.5, 100)


def comparisonPlot(testQuantity, testLabel, controlFile, controlLabel):
    if SHOW_PLOTS:
        set_params(xlabel="$T$",ylabel=testLabel)
        plot_lines(T,testQuantity,xmax=175,marker=None,color=colors[0])
        plot_file(controlFile,style="lines",marker=None,label=controlLabel,color=colors[1])
        plt.show()
        plt.clf()
        clear_legend_labels()


# QM and PDG HRG files
hadrons ,M ,Q ,B ,S ,_ ,g ,w  = readTable("../../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt",
                                           usecols=(0,1,2,3,4,5,6,7),dtype="U11,f8,i8,i8,i8,i8,i8,i8")
hadrons1,M1,Q1,B1,S1,C1,g1,w1 = readTable("../../latqcdtools/physics/HRGtables/PDG_hadron_list_ext_2020.txt",
                                           dtype="U11,f8,i8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))


QMhrg      = HRG(M,g,w,B,S,Q)
pdghrg     = HRG(M1,g1,w1,B1,S1,Q1)
evhrg      = EVHRG(M,g,w,B,S,Q)
evpdghrg   = EVHRG(M1,g1,w1,B1,S1,Q1)
QMhrgexact = HRGexact(M,g,w,B,S,Q)

#
# Test: Calculate chi^200_BQS with b=1 and mu/T=1. Compare against trusted control result.
#
b = 1

chi_QM  = QMhrg.gen_chi(T , B_order=2, Q_order=0, S_order=0, muB_div_T=1)
chi_pdg = pdghrg.gen_chi(T, B_order=2, Q_order=0, S_order=0, muB_div_T=1)
chi_ev  = evhrg.gen_chi(T, b, 1, B_order=2, Q_order=0, S_order=0) + evhrg.gen_chi(T, b, -1, B_order=2, Q_order=0, S_order=0)
chi_ev1 = evpdghrg.gen_chi(T, b, 1, B_order=2, Q_order=0, S_order=0) + evpdghrg.gen_chi(T, b, -1, B_order=2, Q_order=0, S_order=0)

_, refPDG, refQM, refEV, refEV1 = readTable("HRGcontrol/chiBQS_200_muB1.00_b1.00_QMHRG2020_BI.control")

print_results(chi_pdg, refPDG, prec=EPSILON, text="chiB2 PDG check")
print_results(chi_QM , refQM , prec=EPSILON, text="chiB2 QM check")
print_results(chi_ev , refEV , prec=EPSILON, text="chiB2 EV check")
print_results(chi_ev1, refEV1, prec=EPSILON, text="chiB2 EV1 check")

#
# Test : Calculate chi^{101}_uds with b=1 and mu/T=0. Compare against trusted control result.

chiusPDG = pdghrg.genChiFlavor(T,  u_order = 1, d_order = 0, s_order = 1)
chiusQM  = QMhrg.genChiFlavor(T,  u_order = 1, d_order = 0, s_order = 1)


_, refPDGBS, refQMBS = readTable("HRGcontrol/chiBQSC1010_muBh0.0.txt.control")
_, refPDGQS, refQMQS = readTable("HRGcontrol/chiBQSC0110_muBh0.0.txt.control")


print_results(chiusQM / (refQMBS + refQMQS) , np.ones(chiusQM.size)  , prec=EPSILON , text="chius QM check")
print_results(chiusPDG / (refPDGBS + refPDGQS) , np.ones(chiusQM.size) , prec=EPSILON , text="chius PDG check")


#
# Test: Compare some quantities them with PHYSICAL REVIEW D 90, 094503 (2014). We allow for a 2% error
#       tolerance because I had to fish the paper values out by eye and because we are using an updated resonance list
#       compared to what was available in 2014.
#
refT, ref3p_div_T4 = readTable("HRGcontrol/2014_3P_div_T4.d")
test3p_div_T4      = 3*pdghrg.P_div_T4(refT,0,0,0)
print_results(ref3p_div_T4, test3p_div_T4, prec=2e-2, text="2014 HotQCD 3p/T^4 check")
comparisonPlot(3*pdghrg.P_div_T4(T,0,0,0),"$3P/T^4$","HRGcontrol/2014_3P_div_T4.d","2014 HotQCD")


refT, refE_div_T4 = readTable("HRGcontrol/2014_e_div_T4.d")
testE_div_T4      = pdghrg.E_div_T4(refT,0,0,0)
print_results(refE_div_T4, testE_div_T4, prec=3e-2, text="2014 HotQCD e/T^4 check")
comparisonPlot(pdghrg.E_div_T4(T,0,0,0),"$E/T^4$","HRGcontrol/2014_e_div_T4.d","2014 HotQCD")


refT, ref3S_div_4T3 = readTable("HRGcontrol/2014_3s_div_4T3.d")
test3s_div_4T3      = 3*pdghrg.S_div_T3(refT,0,0,0)/4
print_results(ref3S_div_4T3, test3s_div_4T3, prec=3e-2, text="2014 HotQCD 3s/4T^3 check")
comparisonPlot(3*pdghrg.S_div_T3(T,0,0,0)/4,"$3s/4T^3$","HRGcontrol/2014_3s_div_4T3.d","2014 HotQCD")


refT, ref3CV_div_T3 = readTable("HRGcontrol/2014_CV_div_T3.d")
testCV_div_T3       = pdghrg.CV_div_T3_mu0(refT)
print_results(ref3CV_div_T3, testCV_div_T3, prec=3e-2, text="2014 HotQCD CV/T^3 check")
comparisonPlot(pdghrg.CV_div_T3_mu0(T),"$C_V/T^3$","HRGcontrol/2014_CV_div_T3.d","2014 HotQCD")


refT, refcs2 = readTable("HRGcontrol/2014_cs2.d")
cs2 = pdghrg.S_div_T3(refT,0,0,0,0)/pdghrg.CV_div_T3_mu0(refT)
print_results(refcs2, cs2, prec=3e-2, text="2014 HotQCD cs^2 check")


#
# Test: Compare the results from the Taylor series of Bessel functions to the results from the numerical integration.
#       We allow a 30% error tolerance, because for higher temperatures, m/T becomes smaller, making the truncated
#       series less exact. These checks tend to be rather slow, so we parallelize them in a rather naive way.
#
def exactHRGTest(case): # TODO: The integration doesn't seem to be reliable yet...

    refT = np.linspace(40, 170, 10)

    if case == 1:
        testp_div_T4  = pdghrg.P_div_T4(refT,0,0,0)
        exactp_div_T4 = QMhrgexact.P_div_T4(refT,0,0,0)
        print_results(res_true=exactp_div_T4, res=testp_div_T4, prec=1e-1, text="exact p/T^4 check")

    elif case == 2:
        testE_div_T4  = pdghrg.E_div_T4(refT,0,0,0)
        exactE_div_T4 = QMhrgexact.E_div_T4(refT,0,0,0)
        print_results(res_true=exactE_div_T4, res=testE_div_T4, prec=2e-1, text="exact e/T^4 check")

    elif case == 3:
        testNX  = pdghrg.gen_chi(refT,B_order=1,S_order=0,Q_order=0,C_order=0,muB_div_T=1,muQ_div_T=0,muS_div_T=0,muC_div_T=0)
        exactNX = QMhrgexact.number_density(refT,charge='B',muB_div_T=1,muQ_div_T=0,muS_div_T=0,muC_div_T=0)
        print_results(res_true=exactNX, res=testNX, prec=3e-1, text="exact NB")

    elif case == 4:
        testNX  = pdghrg.gen_chi(refT,B_order=0,S_order=1,Q_order=0,C_order=0,muB_div_T=0,muQ_div_T=0,muS_div_T=0.1,muC_div_T=0)
        exactNX = QMhrgexact.number_density(refT,charge='S',muB_div_T=0,muQ_div_T=0,muS_div_T=0.1,muC_div_T=0)
        print_results(res_true=exactNX, res=testNX, prec=3e-1, text="exact NS")

    elif case == 5:
        testNX  = pdghrg.gen_chi(refT,B_order=0,S_order=0,Q_order=1,C_order=0,muB_div_T=0,muQ_div_T=0.2,muS_div_T=0,muC_div_T=0)
        exactNX = QMhrgexact.number_density(refT,charge='Q',muB_div_T=0,muQ_div_T=0.2,muS_div_T=0,muC_div_T=0)
        print_results(res_true=exactNX, res=testNX, prec=1e-1, text="exact NQ")

    elif case == 6:
        tests_div_T3  = pdghrg.S_div_T3(refT,muB_div_T=1,muQ_div_T=0.2,muS_div_T=0.1,muC_div_T=0)
        exacts_div_T3 = QMhrgexact.S_div_T3(refT,muB_div_T=1,muQ_div_T=0.2,muS_div_T=0.1,muC_div_T=0)
        print_results(res_true=exacts_div_T3, res=tests_div_T3, prec=2e-1, text="exact S")

    else:
        pass

#parallel_function_eval(exactHRGTest,[1,2,3,4,5,6],8)

#
# Test: Compare charm results against Physics Letters B 737 (2014) 210–215. I compare with their QMHRG. The tolerance
#       here is even higher. But I expect only to get the right ballpark, since we are not using the same states.
#

data = readTable("../../latqcdtools/physics/HRGtables/hadron_list_ext_strange_charm_2020.txt",
                  usecols=(1,2,3,4,5,6),dtype="f8,i8,i8,i8,i8,i8")

# First exclude all states that have C = 0.
openCharmStates = excludeAtCol(np.array(data),4,0)

# Mesons are those states with B = 0. (Fig. 1 in paper.)
openCharmMesons = restrictAtCol(openCharmStates,2,0)
M, Q, B, S, C, g = openCharmMesons[0], openCharmMesons[1], openCharmMesons[2], openCharmMesons[3], openCharmMesons[4], openCharmMesons[5]
w = np.array([1 if ba==0 else -1 for ba in B])
QMhrg = HRG(M,g,w,B,S,Q,C)
refT, refMesonP_div_T4 = readTable("HRGcontrol/2014_P_Mc.d")
mesonP_div_T4 = QMhrg.P_div_T4(refT,0,0,0)
print_results(mesonP_div_T4, refMesonP_div_T4, prec=1.8e-1, text="2014 HotQCD meson open charm")
comparisonPlot(QMhrg.P_div_T4(T,0,0,0),"$P/T^4$","HRGcontrol/2014_P_Mc.d","2014 open charm meson")

# And the Baryons have B!=0. (Fig. 1 in paper.)
openCharmBaryons = excludeAtCol(openCharmStates,2,0)
M, Q, B, S, C, g = openCharmBaryons[0], openCharmBaryons[1], openCharmBaryons[2], openCharmBaryons[3], openCharmBaryons[4], openCharmBaryons[5]
w = np.array([1 if ba==0 else -1 for ba in B])
QMhrg = HRG(M,g,w,B,S,Q,C)
refT, refBaryonP_div_T4 = readTable("HRGcontrol/2014_P_Bc.d")
baryonP_div_T4 = QMhrg.P_div_T4(refT,0,0,0)
print_results(baryonP_div_T4, refBaryonP_div_T4, prec=2.2e-1, text="2014 HotQCD baryon open charm")
comparisonPlot(QMhrg.P_div_T4(T,0,0,0),"$P/T^4$","HRGcontrol/2014_P_Bc.d","2014 open charm baryon")

# Do another check against Frithjof's data.
refT, _, _, refBaryonP_div_T4 = readTable("HRGcontrol/OUT_5.0.DAT140_2022_hidden_charm_pressure")
refT *= 1000 # He gives his temperatures in [GeV]
baryonP_div_T4 = QMhrg.P_div_T4(refT,0,0,0)
print_results(baryonP_div_T4, refBaryonP_div_T4, prec=1e-3, text="2022 F. Karsch code open charm")


# Finally we check one of the derivatives. (Fig. 4 in paper.)
M, Q, B, S, C, g = openCharmStates[0], openCharmStates[1], openCharmStates[2], openCharmStates[3], openCharmStates[4], openCharmStates[5]

w = np.array([1 if ba==0 else -1 for ba in B])
QMhrg = HRG(M,g,w,B,S,Q,C)
refT, refRSC13 = readTable("HRGcontrol/2014_chi_charm_ratio.d")
chiBSC112 = QMhrg.gen_chi(refT,B_order=1,S_order=1,Q_order=0,C_order=2)
chiSC13   = QMhrg.gen_chi(refT,B_order=0,S_order=1,Q_order=0,C_order=3)
RSC13     = -chiBSC112/(chiSC13 - chiBSC112)
print_results(RSC13, refRSC13, prec=1.4e-1, text="2014 HotQCD RSC13")


#
# Test: Compare numerical derivatives against analytic derivatives.
#
T   = np.linspace(10,150,140)

muBList = [1,3,10,30,100]
muSList = [1,3,10,30,100]
for muBh in muBList:
    for muSh in muBList:

        logger.info("TESTS AT muB/T, muS/T = ",muBh,muSh)

        left  = QMhrg.P_div_T4(T,muB_div_T=muBh,muS_div_T=muSh)
        right = QMhrg.gen_chi(T,B_order=0,Q_order=0,S_order=0,muB_div_T=muBh,muS_div_T=muSh)
        print_results(left, right, prec=EPSILON, text="P vs gen_chi")

        def Ehat(t):
            return QMhrg.E_div_T4(t,muB_div_T=muBh,muS_div_T=muSh)
        exact     = QMhrg.ddT_E_div_T4(T, muB_div_T=muBh,muS_div_T=muSh)
        numerical = diff_deriv(T,Ehat)
        print_results(exact, numerical, prec=EPSILON, text="d(E/T^4)/dT")

        def Phat(t):
            return QMhrg.P_div_T4(t,muB_div_T=muBh,muS_div_T=muSh)
        exact     = QMhrg.ddT_P_div_T4(T, muB_div_T=muBh,muS_div_T=muSh)
        numerical = diff_deriv(T,Phat)
        print_results(exact, numerical, prec=EPSILON, text="d(P/T^4)/dT")

        def Shat(t):
            return QMhrg.S_div_T3(t,muB_div_T=muBh,muS_div_T=muSh)
        exact     = QMhrg.ddT_S_div_T3(T, muB_div_T=muBh,muS_div_T=muSh)
        numerical = diff_deriv(T,Shat)
        print_results(exact, numerical, prec=3*EPSILON, text="d(S/T^4)/dT")

        def chiBQ11(t):
            return QMhrg.gen_chi(t,B_order=1,Q_order=1,S_order=0,C_order=0,muB_div_T=muBh,muS_div_T=muSh)
        exact     = QMhrg.ddT_gen_chi(T,B_order=1,Q_order=1,S_order=0,C_order=0,muB_div_T=muBh,muS_div_T=muSh)
        numerical = diff_deriv(T,chiBQ11)
        print_results(exact, numerical, prec=EPSILON, text="d(chi11BQ)/dT")

        def ddT_chiBQ11(t):
            return QMhrg.ddT_gen_chi(t, B_order=1, Q_order=1, S_order=0, C_order=0, muB_div_T=muBh,muS_div_T=muSh)
        exact     = QMhrg.d2dT2_gen_chi(T,B_order=1,Q_order=1,S_order=0,C_order=0, muB_div_T=muBh,muS_div_T=muSh)
        numerical = diff_deriv(T,ddT_chiBQ11)
        print_results(exact, numerical, prec=EPSILON, text="d^2(chi11BQ)/dT^2")


muh = np.linspace(0,1.5,len(T))


exact = QMhrg.gen_ddmuh_E_div_T4(T,B_order=1,Q_order=0,S_order=0,C_order=0,muB_div_T=muh)
def Ehat(muh):
    return QMhrg.E_div_T4(T,muB_div_T=muh)
numerical = diff_deriv(muh,Ehat)
print_results(exact, numerical, prec=1e-4, text="d(E/T^4)/dmuB")


exact = QMhrg.gen_ddmuh_P_div_T4(T,B_order=1,Q_order=0,S_order=0,C_order=0,muB_div_T=muh)
def Phat(muh):
    return QMhrg.P_div_T4(T,muB_div_T=muh)
numerical = diff_deriv(muh,Phat)
print_results(exact, numerical, prec=1e-4, text="d(P/T^4)/dmuB")


exact = QMhrg.gen_ddmuh_S_div_T3(T,B_order=1,Q_order=0,S_order=0,C_order=0,muB_div_T=muh)
def Shat(muh):
    return QMhrg.S_div_T3(T,muB_div_T=muh)
numerical = diff_deriv(muh,Shat)
print_results(exact, numerical, prec=1e-4, text="d(S/T^4)/dmuB")


times.printTiming()
