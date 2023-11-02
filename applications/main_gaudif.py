#!/bin/python3

# 
# main_gaudif.py                                                               
# 
# D. Clarke 
# 
# Do a Gaussian difference test on user-provided data.
# 

import argparse
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import gaudif
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.utilities import getArgs

parser = argparse.ArgumentParser(description='Do a Z-test.')
parser.add_argument('--xm1', dest='xm1', required=True, type=float, help='mean 1')
parser.add_argument('--eb1', dest='eb1', required=True, type=float, help='error bar 1')
parser.add_argument('--xm2', dest='xm2', required=True, type=float, help='mean 2')
parser.add_argument('--eb2', dest='eb2', required=True, type=float, help='error bar 2')

args=getArgs(parser)

xm1=args.xm1
eb1=args.eb1
xm2=args.xm2
eb2=args.eb2

q=gaudif(xm1,eb1,xm2,eb2)

logger.info(" meas1 = ",get_err_str(xm1,eb1))
logger.info(" meas2 = ",get_err_str(xm2,eb2))
logger.info("     q = ",q)
