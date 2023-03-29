import numpy as np
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger
from latqcdtools.physics.HotQCDEOS import EOS
from latqcdtools.base.readWrite import readTable, writeTable


Nnparams = readTable("../latqcdtools/physics/EoSparams/Eos_params_Nn2022.txt", unpack=True)
qnparams = readTable("../latqcdtools/physics/EoSparams/Eos_params_qn2022.txt", unpack=True)

T = np.arange(135,285,5)
eos = EOS(T)

muBdivT = 2.0

p, nB, e, s = eos.ObsEoSfixedmuB(Nnparams, qnparams, muBdivT)

outFileName = "HotQCDEos2022_fixedmuBdivT%0.2f.txt"%muBdivT
writeTable(outFileName, T, muBdivT * T , p, nB, e , s, header='TMeV muBMeV p nB e s')


snB = 100

muBdivT, p, nB, e, s = eos.ObsEoSfixedsnB(Nnparams, qnparams, snB)

outFileName = "HotQCDEos2022_fixedsnB%0.2f.txt"%snB
writeTable(outFileName, T, muBdivT * T , p, nB, e , s, header='TMeV muBMeV p nB e s')
