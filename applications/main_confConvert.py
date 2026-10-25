#!/bin/python3

#
# main_confConvert.py
#
# D. Clarke
#
# Convert a gauge configuration from one format to another. The format and lattice size of the input
# configuration are detected from the file.
#

import argparse
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import getArgs
from latqcdtools.interfaces.confReader import detectConf, NERSCReader, ILDGReader, DynisoReader
from latqcdtools.interfaces.confWriter import NERSCWriter

logger.set_log_level('INFO')

# To support another output format, implement its writer in confWriter.py and add it here.
readers = { 'NERSC'  : NERSCReader,
            'ILDG'   : ILDGReader,
            'dyniso' : DynisoReader }
writers = { 'NERSC'  : NERSCWriter }

parser = argparse.ArgumentParser(description='Convert a gauge configuration to another format.')
parser.add_argument('--inconf', dest='inconf', required=True, help='input configuration (NERSC, ILDG, or dyniso)')
parser.add_argument('--outconf', dest='outconf', required=True, help='output configuration')
parser.add_argument('--outformat', dest='outformat', required=True, type=str.upper, choices=list(writers.keys()),
                    help='format of the output configuration')

args = getArgs(parser)

inconf    = args.inconf
outconf   = args.outconf
outformat = args.outformat

informat, Ns, Nt = detectConf(inconf)
logger.info('Detected',informat,'configuration of size',str(Ns)+'^3 x',Nt)

reader = readers[informat](Ns=Ns, Nt=Nt)
gauge  = reader.readConf(inconf)

writer = writers[outformat](gauge)
writer.writeConf(outconf)
logger.info('Wrote',outformat,'configuration',outconf)
