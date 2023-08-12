#!/bin/python3

#
# J. Goswami
#
# Evaluate EOS from the parametrization given in 2212.10016 [hep-lat] and 1407.6387 [hep-lat]
# The Eos parametrization is currently reliable in, T = [135:300] and muB / T = [0:2.5]
#
# You can use this code by, for instance, 
#
#   python main_HotQCDEoS.py --EosType "fixedsnB" --snB 50
#   python main_HotQCDEoS.py --EosType "fixedmuB" --muBdivT 2.0
#
# Usage help: python main_HotQCDEoS.py -h
#

from latqcdtools.physics.HotQCDEOS import EOS
import numpy as np
import argparse
from latqcdtools.base.utilities import getArgs
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable, writeTable

logger.set_log_level('INFO')

parser = argparse.ArgumentParser(description='Script to calculate '
                                             'pressure, energy density and entropy density '
                                             'from the parametrization given in 2212.10016 [hep-lat] for r = 0.4', allow_abbrev=False)
parser.add_argument("--EosType", dest="EoT", help="choose one:'fixedmuB' or 'fixednB'", required=True, type=str)
parser.add_argument("--muBdivT", dest="muBdivT", default=None, help="If 'fixedmuB' option choosen : Value of constant muBdivT", type=float)
parser.add_argument("--snB", dest="snB", default=None, help="If 'fixedsnB' option choosen : Value of constant s/nB", type=float)

# parametrization of the number density and q's
Nnparams = readTable("../latqcdtools/physics/EoSparams/Eos_params_Nn2022.txt", unpack=True)
qnparams = readTable("../latqcdtools/physics/EoSparams/Eos_params_qn2022.txt", unpack=True)

T = np.arange(135, 285, 5)
eos = EOS(T)

args = getArgs(parser)

muBdivT = args.muBdivT
snB = args.snB
eot = args.EoT

if muBdivT is not None and eot == "fixedmuB":
    p, nB, e, s = eos.ObsEoSfixedmuB(Nnparams, qnparams, muBdivT)
    outFileName = "HotQCDEos2022_fixedmuBdivT%0.2f.txt" % muBdivT
    writeTable(outFileName, T, muBdivT * T, p, nB, e, s, header='TMeV muBMeV p nB e s')
elif snB is not None and eot == "fixedsnB":
    muBdivT, p, nB, e, s = eos.ObsEoSfixedsnB(Nnparams, qnparams, snB)
    outFileName = "HotQCDEos2022_fixedsnB%0.2f.txt" % snB
    writeTable(outFileName, T, muBdivT * T, p, nB, e, s, header='TMeV muBMeV p nB e s')
elif muBdivT == 0.0:
    p, e, s = eos.ObsEoS()
    outFileName = "HotQCDEos2022_fixedmuBdivT%0.1f.txt" % muBdivT
    writeTable(outFileName, T, muBdivT * T, p, e, s, header='TMeV muBMeV p nB e s')
else:
    print("To calculate EOS either provide a fixed muB/T value or fixed s/nB value see : python  main_HotQCDEoS.py -h")
