#!/bin/bash

BQS1=200
BQS2=130
filepath="../latqcdtools/physics/QM_hadron_list_ext_strange_2020.txt" # change this line accordingly !!
b=1.0 #excluded volume parameter

#Here we demonstrate two usage of the hrg 

# col number should provided in the following way
# hadron  = Name of the hadron [All particles and anti-particles]
# M       = Mass of the hadron [in MeV recommended]
# Q       = charge of the hadron
# B       = Baryon number of the hadron [B=1,-1,0,0 for baryon,anti-baryon,meson,anti-mesons]
# S       = Strangeness number of the hadron
# C       = Charm number of the hadron
# g       = degenracy of the hadron state

python3 main_evhrg.py --hadron_file ${filepath} --bqs ${BQS1} --b $b --column 0 1 2 3 4 5 6 

# temperature in MeV. start:end:interval
python3 main_evhrg.py --hadron_file ${filepath} --bqs ${BQS2} --b $b --column 0 1 2 3 4 5 6 --temperature_range 130:180:0.5
