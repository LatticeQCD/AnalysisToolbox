r=0.4 #nQ/nB=0.4

filepath="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt"

time python3 main_HRG_LCP.py --r $r --hadron_file ${filepath} --models QM --T 100
