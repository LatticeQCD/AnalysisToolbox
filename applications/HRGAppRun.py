#!/bin/python3

# 
# J. Goswami 
# 
# This is an auxiliary code that shows how you can use the applications
# main_HRG_LCP or main_HRG_measure to do production runs. 
# 

import argparse
import latqcdtools.base.logger as logger
from latqcdtools.base.initialize import initialize, finalize
from latqcdtools.base.speedify import parallel_function_eval
from latqcdtools.base.utilities import getArgs, shellVerbose


def task(temp, mode):
    """ Speed things up with parallelization. """
    BQSC = "2000"
    if mode == 1:
        shellVerbose("python3", "main_HRG_LCP.py", "--r", r, "--hadron_file", filepath, "--models", "QM", "--T", temp)
    elif mode == 2:
        shellVerbose("python3", "main_HRG_measure.py", "--hadron_file", filepath,
                     "--LCP_file", "HRG_LCP_T%0.1f_r0.4QM"%temp, "--bqsc", BQSC, "--obs", "chi")


parser = argparse.ArgumentParser(description='Run different modes of the script.')
parser.add_argument('runMode', type=int, choices=[1, 2, 3], help='The run mode of the script.')
args = getArgs(parser)

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
                shellVerbose("python3", "main_strangeness_neutral.py", "--hadron_file", filepath, "--fixedmuBNszerofile",
                             f"fixedmuBTNsZerofiles/HRG_fixedmuBT{mb}_r0.4QMHRG2020_BI", "--obs", obs, "--r", r, "--tag", "QMHRG2020_BI_Nszero")
    else:
        logger.TBRaise("Invalid runMode",args.runMode)

    finalize()

if __name__ == '__main__':
    main()
