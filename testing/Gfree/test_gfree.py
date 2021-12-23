#!/usr/bin/env python3
from latqcdtools.gfree import *
from latqcdtools.tools import *
import numpy as np


nts, Gs = np.loadtxt("ref.txt").transpose()

nt = 15
mdT = 14
#From scalar
G_ref = Gfree_vec(15/len(nts), mdT)
print("nt = ", nt, "Agreement with reference value?",
        Gs[nt], rel_check(Gs[nt], G_ref))
#From array
Gs_ref = Gfree_vec(nts/len(nts), mdT)


for nt in range(1, len(nts)):
    print("nt = ", nt, "Agreement with reference value?",
            Gs[nt], rel_check(Gs[nt], Gs_ref[nt]))

