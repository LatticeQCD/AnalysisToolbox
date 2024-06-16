latqcdtools.physics.HRG
=============

`LCP_init_NS0(muB)`
 
    Give a good initial guess for NS=0 LCP. 
    
`RMS_mass(Nt, T)`


`dmuh(order, muh)`
 
    d^order/dmuh^order derivative of muh 
    
`EVHRG(Mass, g, w, B, S, Q, C=None)`
 
    Excluded volume hadron resonance gas. b is excluded volume parameter. 
    
`HRG(Mass, g, w, B, S, Q, C=None, NMAX_light=21, NMAX_heavy=2)`
 
    HRG implemented through Taylor expasion of logarithm. For more information please see, e.g.
    Physics Letters B 695 (2011) 136–142 or especially arXiv:2011.02812.
    You can optionally adjust NMAX_light and NMAX_heavy, which control the number of terms to keep in the
    Taylor expansion for species that are respectively lighter and heavier than the Kaon.
    Our pressure is given in terms of a Taylor series involving modified Bessel functions of the second kind, which
    needs to be truncated at some order. These functions get strongly suppressed when their argument is large.
    In our case, this argument is proportional to the mass. Hence we will sometimes take the Boltzmann
    approximation, i.e. that the mass is large compared to the temperature. In this limit, even fewer terms of the
    expansion need to be kept. Doing so boosts performance.

    We work in the grand canonical ensemble, i.e. we consider P(V,T,mu_i/T). Hence in this class, derivatives w.r.t.
    one of those N_chemical_potentials + 2 variables assume all others are held fixed. 
    
`HRGbase(Mass, g, w, B, S, Q, C=None)`
 
    Hadron resonance gas base class. Here we collect methods and attributes that all HRG-type classes should have
    in common. Mass=mass of the Hadron/resonance , g=spin degenerecy , w= fermi(-1)/bose(1) statistics.
    B, Q, S, and C are respectively the baryon number, electric charge, strangeness, and charm of each state. 
    
`HRGexact(Mass, g, w, B, S, Q, C=None)`
 
    HRG implemented through numerical integration. 
    
