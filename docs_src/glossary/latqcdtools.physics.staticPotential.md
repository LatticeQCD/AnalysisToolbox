latqcdtools.physics.staticPotential
=============

`V_Teq0(r) -> float`
 
    Zero temperature quark potential in [MeV], takes r in [fm]. The parameters a, b, and c come from
    a Levenberg-Marquardt fit of the data in Fig 14 of Phys. Rev. D90 (2014) 094503. These numbers can
    be obtained again by running analysistoolbox/hisq_potential/fit_hisq_pot.py. 
    
`fitV_Teq0(r, a, b, c) -> float`
 
    Fit form of standard Cornell potential. Fit to be done in lattice units.
    
`fitV_Teq0_oneloop(r, a, b, c, d) -> float`
 
    Including one-loop corrections to Coulomb. See Nucl. Phys. B 129 (1977) and Phys. Lett B 92 (1980).
    Fit to be done in lattice units.
    
`fitV_Teq0_twoloop(r, a, b, c, d, e) -> float`
 
    Including two-loop corrections to Coulomb. See Nucl. Phys. B 501 (1997). Fit to be done in lattice units.
    
`impdist(Ns, r2max, improvedAction=True)`

    Calculation of tree-level improved distances. Follows eq. (3) of 10.1103/PhysRevD.90.074038,

    INPUT:
           Ns--spatial extension of lattice
        r2max--maximum squared distance to improve

    OUTPUT:
        rimp--list of improved distances
    
