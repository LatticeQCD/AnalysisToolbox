#!/bin/bash


filepath="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt"

filepath_charm="../latqcdtools/physics/HRGtables/hadron_list_ext_strange_charm_2020.txt"

# hadron  = Name of the hadron [All particles and anti-particles]
# M       = Mass of the hadron [in MeV recommended]
# Q       = charge of the hadron
# B       = Baryon number of the hadron [B=1,-1,0,0 for baryon,anti-baryon,meson,anti-mesons]
# S       = Strangeness number of the hadron
# C       = Charm number of the hadron
# g       = degenracy of the hadron state


# This part is for examples/testing
BQSC=2000
python3 main_HRG_measure.py --hadron_file ${filepath} --tag QMHRG2020_BI --obs chi --bqsc ${BQSC} --temperature_range 130:180:0.5 --muB 0.0

