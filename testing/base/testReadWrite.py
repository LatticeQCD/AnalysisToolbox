# 
# testReadWrite.py                                                               
# 
# D. Clarke
# 
# Tests for read and write methods.
#
import argparse
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import writeTable
from latqcdtools.base.utilities import getArgs

parser = argparse.ArgumentParser(description='Tests for read and write methods.')

parser.add_argument('--type', dest='headerType', type=str, default='str', help='Is header a list or string?')

args = getArgs(parser)

xdata = np.linspace(0,1,101)
ydata = xdata**2

if args.headerType=="str":
    writeTable("test.d",xdata,ydata,header="         who         bastank")
elif args.headerType=="list":
    writeTable("test.d",xdata,ydata,header=["who","bastank"])
else:
    logger.TBError("Invalid header type",args.headerType)