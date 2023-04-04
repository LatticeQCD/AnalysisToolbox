import numpy as np
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger
from latqcdtools.physics.HotQCDEOS import EOS
from latqcdtools.base.readWrite import readTable


logger.set_log_level('INFO')

EPSILON = 1e-6
refT, refp, refe, refs, _ = readTable("EoScontrol/Thermodynamic_observables_zero_chemical_potential_T130-300.txt",
                                      unpack=True)
Nnparams = readTable("../../latqcdtools/physics/EoSparams/Eos_params_Nn2022.txt", unpack=True)
qnparams = readTable("../../latqcdtools/physics/EoSparams/Eos_params_qn2022.txt", unpack=True)

T = refT
eos = EOS(T)

# pressure , energy density and entropy density test for muB = 0, ref: arXiv : 1407.6387
# If anyone wants the muB  dependent part only one needs to subtract the p0,e0 and s0
p0, e0, s0 = eos.ObsEoS()

muB = 0.0
p, nB, e, s = eos.ObsEoSfixedmuB(Nnparams, qnparams, muB)

print_results(p / p0, np.ones(refT.size), prec=EPSILON / 1e6,
              text="Self consistency check")
print_results(e / e0, np.ones(refT.size), prec=EPSILON / 1e6,
              text="Self consistency check")

print_results(s / s0, np.ones(refT.size), prec=EPSILON / 1e6,
              text="Self consistency check")

print_results(p / refp, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the pressure for muB = 0 Test (ref : arXiv : 1407.6387)")
print_results(e / refe, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the energy desnity for muB = 0 Test (ref : arXiv : 1407.6387)")
print_results(s / refs, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the entropy desnity for muB = 0 Test (ref : arXiv : 1407.6387)")

# pressure , energy density and entropy density test for muB = 1.0, ref: arXiv:2212.10016 [hep-lat], 2212.09043 [
# hep-lat], 1701.04325 [hep-lat]

EPSILON = 1e-5
muB = 1.0  # Also try muB = 1.0

p, nB, e, s = eos.ObsEoSfixedmuB(Nnparams, qnparams, muB)

refT, refp, refe, refs = readTable("EoScontrol/Eos_control_muB%0.1f.txt" % muB, unpack=True)

print_results(p / refp, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the pressure for muB = %0.1f Test (ref : arXiv : 2212.10016)" % muB)
print_results(e / refe, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the energy desnity for muB = %0.1f Test (ref : arXiv : 2212.10016)" % muB)
print_results(s / refs, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the entropy desnity for muB = %0.1f Test (ref : arXiv : 2212.10016)" % muB)

# pressure , energy density and entropy density test for s/nB = 50 , ref: arXiv:2212.10016 [hep-lat], 2212.09043 [
# hep-lat]

snB = 50  # Also try snB = 200

muB, p, nB, e, s = eos.ObsEoSfixedsnB(Nnparams, qnparams, snB)

refT, refmuB, refp, refe, refs = readTable("EoScontrol/Eos_control_snB%d.txt" % snB, unpack=True)

print_results(muB / refmuB, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the muB trajectory for snB = %d Test (ref : arXiv : 2212.10016)" % snB)
print_results(p / refp, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the pressure for snB = %d Test (ref : arXiv : 2212.10016)" % snB)
print_results(e / refe, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the energy desnity for muB = %d Test (ref : arXiv : 2212.10016)" % snB)
print_results(s / refs, np.ones(refT.size), prec=EPSILON,
              text="Given parametrization of the entropy desnity for muB = %d Test (ref : arXiv : 2212.10016)" % snB)
