#!/bin/bash
b=1.0 # value of the excluded volume parameter !!
Tpc=156.5 # Tpc value in MeV
r=0.4 #nQ/nB=0.4


filepath="../latqcdtools/physics/QM_hadron_list_ext_strange_2020.txt"
filepath_muB="Tpc156.5_pseudo-muS_muB_r0.4QMHRG2020_BI"

#python3 strangeness_neutral_hrg.py --Tpc ${Tpc} --r $r --hadron_file ${filepath} --b $b --tag QMHRG2020_BI
python3 pseudo_sus_evhrg.py --hadron_file ${filepath} --tag QMHRG2020_BI --bqs 100 --mubmus_file ${filepath_muB} --b ${b} --r 0.4

# This is exact solution of muS and muB

