#!/bin/bash


r=0.4 #nQ/nB=0.4
filepath="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt"


# 0: Generate LCP at many temperatures.
# 1: Make measurements along LCP.
# 2: Measure specifically cs2.
runMode=2


if [ ${runMode} -eq 0 ]; then

  temps=($(seq 30 1 160))
  for temp in "${temps[@]}"; do
    python3 main_HRG_LCP.py --r $r --hadron_file ${filepath} --models QM --T ${temp}
  done

elif [ ${runMode} -eq 1 ]; then

  python3 main_HRG_measure.py --hadron_file ${filepath} --models QM --LCP_file HRG_LCP_T100.0_r0.4 --bqsc 2000 --obs cs2

elif [ ${runMode} -eq 2 ]; then

  python3 main_HRG_cs2.py --hadron_file ${filepath}  --LCP_file HRG_LCP_T100.0_r0.4 

else

  echo "Invalid runMode"
  exit

fi
