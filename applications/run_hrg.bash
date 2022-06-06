#!/bin/bash


filepath="../latqcdtools/physics/QM_hadron_list_ext_strange_2020.txt" # change this line accordingly !!


# col number should provided in the following way
# hadron  = Name of the hadron [All particles and anti-particles]
# M       = Mass of the hadron [in MeV recommended]
# Q       = charge of the hadron
# B       = Baryon number of the hadron [B=1,-1,0,0 for baryon,anti-baryon,meson,anti-mesons]
# S       = Strangeness number of the hadron
# C       = Charm number of the hadron
# g       = degenracy of the hadron state


# This part is for examples/testing
BQS1=200
BQS2=200
b=1 #excluded volume parameter
#python3 main_evhrg.py --hadron_file ${filepath} --tag QMHRG2020_BI --obs chi --bqs ${BQS1} --b $b  
# temperature in MeV. start:end:interval
python3 main_evhrg.py --hadron_file ${filepath} --tag QMHRG2020_BI --obs chi --muB 1.0 --bqs ${BQS2} --b $b --temperature_range 130:180:0.5


# Here's where you can actually use it.
#python3 main_evhrg.py --hadron_file ../latqcdtools/physics/QM_hadron_list_ext_strange_2020.txt --tag QMHRG2020_BI --obs chi --bqs 101 --b 0.4 --temperature_range 130:180:0.5
#python3 main_muBhrg.py --hadron_file ../latqcdtools/physics/QM_hadron_list_ext_strange_2020.txt --tag QMHRG2020_BI --obs chi --bqs 101 --b 0.4 --temperature_range 130:180:0.5
