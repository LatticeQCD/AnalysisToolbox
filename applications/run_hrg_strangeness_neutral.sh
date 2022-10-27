#!/bin/bash


r=0.5 #nQ/nB
filepath="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt"


# 0: Generate LCP at many temperatures.
# 1: Measure observables. 
# 2: Measure observables at fixed muB/T with Ns=0 .
runMode=0
NTASKS=4

if [ ${runMode} -eq 0 ]; then

  temps=($(seq 16 16))

  task() {
    python3 main_HRG_LCP.py --r $r --hadron_file ${filepath} --models QM --T $1
  }

  for temp in "${temps[@]}"; do
    ((i=i%NTASKS)); ((i++==0)) && wait
    task "${temp}" &
  done

elif [ ${runMode} -eq 1 ]; then

  python3 main_HRG_measure.py --hadron_file ${filepath} --models QM --LCP_file HRG_LCP_T100.0_r0.4 --bqsc 2000 --obs chi 

elif [ ${runMode} -eq 2 ]; then

  for obs in energy specificheat; do
    for mb in 0.0 1.0 1.5 2.0 2.5; do
      python3 main_strangeness_neutral.py --hadron_file ${filepath}  --fixedmuBNszerofile HRG_fixedmuBT${mb}_r0.4QMHRG2020_BI --obs ${obs} --r $r --tag QMHRG2020_BI_Nszero 
    done
  done

else
  echo "Invalid runMode"
  exit
fi
