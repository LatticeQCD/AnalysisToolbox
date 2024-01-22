#!/bin/python3

# 
# J. Goswami 
# 
# This is an auxiliary code that shows how you can use the applications
# main_HRG_LCP or main_HRG_measure to do production runs. 
# 

import subprocess
import argparse
from latqcdtools.base.initialize import initialize, finalize
from latqcdtools.base.speedify import parallel_function_eval


def task(temp, mode):
    """ Speed things up with parallelization. """
    BQSC = "2000"
    if mode == 1:
        subprocess.run(["python3", "main_HRG_LCP.py", "--r", str(r),
                       "--hadron_file", filepath, "--models", "QM", "--T", str(temp)])
    elif mode == 2:
        subprocess.run(["python3", "main_HRG_measure.py", "--hadron_file", filepath,
                       "--LCP_file", "HRG_LCP_T%0.1f_r0.4QM"%temp, "--bqsc", BQSC, "--obs", "chi"])


parser = argparse.ArgumentParser(description='Run different modes of the script.')
parser.add_argument('runMode', type=int, choices=[1, 2, 3], help='The run mode of the script.')
args = parser.parse_args()

r = 0.4
filepath = "../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt"
temps = range(120, 135, 5)


def main():

    initialize('HRG.log')

    if args.runMode == 1:
        parallel_function_eval(task, temps, args=([1])) 

    elif args.runMode == 2:
        parallel_function_eval(task, temps, args=([2])) 

    elif args.runMode == 3:
        for obs in ["pressure", "energy"]:
            for mb in [0.0, 1.0, 1.5, 2.0, 2.5]:
                subprocess.run(["python3", "main_strangeness_neutral.py", "--hadron_file", filepath, "--fixedmuBNszerofile",
                               f"fixedmuBTNsZerofiles/HRG_fixedmuBT{mb}_r0.4QMHRG2020_BI", "--obs", obs, "--r", str(r), "--tag", "QMHRG2020_BI_Nszero"])
    else:
        print("Invalid runMode")

    finalize()

if __name__ == '__main__':
    main()
