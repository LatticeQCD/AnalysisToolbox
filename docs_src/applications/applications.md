Applications
=============

We have a few ready-to-use applications for some basic physics calculations. If you make a new application,
please make sure it starts with `main_`.

## main_getBeta.py
Given a reference scale, $N_\tau$, and $T$ in [MeV], this script calculates the corresponding $\beta$. An
example calling is
```shell
./main_getTempandSpacing.py --Nt 8 --scale r0 --T 300 
```

## main_getTempAndSpacing.py
Given a reference scale, $N_\tau$, and $\beta$ value, this script calculates $T$ in [MeV] and $a$ in [fm]. An
example calling is
```shell
./main_getTempandSpacing.py --Nt 8 --beta 6.285 --scale fk
```
For more details about the scale setting, please look into the 
[reference scales](../physicsAnalysis/referenceScales.md) article.

## main_HRG_measure.py
The goal of this script is to calculate some observables in the [HRG model](../physicsAnalysis/HRG.md) given some
external control parameters, like a range of temperatures with fixed $\mu_B/T$. By default, this program uses
a list of hadrons and resonances created by the HotQCD collaboration called QMHRG2020. You can find this list here:
```shell
latqcdtools/physics/HRGtables/hadron_list_ext_strange_2020.txt
```
A list that includes charmed hadrons and resonances is
```shell
latqcdtools/physics/HRGtables/hadron_list_ext_strange_charm_2020.txt
```
You can choose your hadron list with the `--hadron_file` argument. It is also up to you to choose a model; at the moment
the only possibilities are the typical HRG `QM` or an excluded volume HRG model `EV`. This is specified with
the argument `--model`. Finally you need to specify your observable `--obs`. If you specify a generalized
susceptibility `chi`, you must also pass the $B$, $Q$, $S$, and $C$ chemical potential derivative orders.
A straightforward usage of this script is, for instance,
```shell
./main_HRG_measure.py --obs chi --bqsc 1100 --temperature_range 130:180:0.5 --muB 0.0
```

## main_HRG_LCP.py
In the context of QCD at finite chemical potential, it is interesting to examine systems following a few lines
of constant physics (LCPs). The HotQCD collaboration has focused on strangeness-neutral systems, with $n_Q/n_B=0.4$
(corresponding to gold-gold collisions at RHIC) and $n_Q/n_B=0.5$ (corresponding to the isospin-symmetric case).
This script creates tables of $\mu_B/T$, $\mu_Q/T$, and $\mu_S/T$ ($\mu_C=0$) that lie on $n_S=0$ LCPs like the ones
mentioned above.
A straightforward usage of this script is, for instance,
```shell
./main_HRG_LCP.py --r 0.4 --models QM --T 150
```
Once you have generated some LCP files, you can also use `main_HRG_measure.py` from above to carry out measurements
on them. In such a case, you must pass the LCP file as argument. For instance,
```shell
./main_HRG_measure.py --obs chi --bqsc 1100 --LCP_file HRG_LCP_T150.0_r0.5QM
```

## main_plotRatApprox.py
The RHMC of [SIMULATeQCD](https://github.com/LatticeQCD/SIMULATeQCD) relies on a rational approximation to the fermion
determinant. It is useful to see how well this approximation compares with the exact function. This can be checked
visually with, for example
```shell
./main_plotRatApprox.py in.rational 0.001 0.01
```
where `in.rational` is the rational approximation file, and we use a light quark mass of 0.001 in lattice units
and a strange quark mass of 0.01.

## main_HotQCDEoS.py 
The paramterization for HotQCD equation of state at  $\mu_B/T = 0$, $\mu_Q/T = 0$, and $\mu_S/T = 0$ are given in "Equation of state in ( 2+1 )-flavor QCD, Phys.Rev.D 90 (2014) 094503, (HotQCD Collaboration) A. Bazavov et al.". The pressure ($P$), energy density ($\epsilon$) and entropy density ($s$) can be obtained from the thermodynamic relations.
```shell
./main_HotQCDEoS.py --EosType "fixedmuB" --muBdivT 0.0
```
