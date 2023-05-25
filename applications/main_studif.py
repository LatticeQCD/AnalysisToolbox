#!/bin/python3

# 
# main_studif.py                                                               
# 
# D. Clarke 
# 
# Do a Student difference test on user-provided data.
# 

import argparse
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import studif
from latqcdtools.base.printErrorBars import get_err_str

parser = argparse.ArgumentParser(description='Do a Z-test.')
parser.add_argument('--xm1', dest='xm1', required=True, type=float, help='mean 1')
parser.add_argument('--eb1', dest='eb1', required=True, type=float, help='error bar 1')
parser.add_argument('--n1' , dest='n1' , required=True, type=int  , help='number of measurements 1')
parser.add_argument('--xm2', dest='xm2', required=True, type=float, help='mean 2')
parser.add_argument('--eb2', dest='eb2', required=True, type=float, help='error bar 2')
parser.add_argument('--n2' , dest='n2' , required=True, type=int  , help='number of measurements 2')

args=parser.parse_args()

xm1=args.xm1
eb1=args.eb1
xm2=args.xm2
eb2=args.eb2
n1 =args.n1
n2 =args.n2

q=studif(xm1,eb1,n1,xm2,eb2,n2)

logger.info(" meas1 = ",get_err_str(xm1,eb1))
logger.info(" meas2 = ",get_err_str(xm2,eb2))
logger.info("     q = ",q)
