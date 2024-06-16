latqcdtools.interfaces.HotQCD
=============

`_badNf(Nf)`


`_badNt(Nt, Nf)`


`_badbeta(beta, msml, Nt, Nf)`


`_badmsml(msml, Nt, Nf)`


`loadDens(densFile, confID, lp, inTable=None) -> dict`
 
    Allows reading of output from C. Schmidt's Dense code. The Dense code produces measurements of various operators
    relevant calculating conserved charge fluctuations. We store as a dictionary indexed by confID, which lets us
    conveniently combine data when there are multiple dense files per configuration. Here we update the table
    inTable for the new read-in, which yields outTable. Only supports Nf=2+1 configurations for the time being.

    Parameters
    ----------
    densFile : str
        Name of the Dense code file to be read.
    confID : str
        A unique identifier for the configuration on which dense calculation was performed. Every person seems to have a
        different scheme for labelling configurations, so this needs to be a string to be as flexible as possible.
    lp : latticeParams
        Parameters for the ensemle the configuration belongs to.
    inTable : dict
        A table indexed by confID. Its values are a list of operators that have been measured.

    Returns
    -------
    outTable : dict
        A table indexed by confID. Its values are a list of operators that have been measured.
    
`makeConfTag(conf, stream) -> str`
 
    This takes a configuration number conf and stream label stream to make a tag labelling a configuration.
    Implementing this as a function makes sure everyone using the Toolbox as the same convention and, more importantly,
    ensures that the tags have no whitespace in them, which otherwise can throw off the column counting of methods in
    the denseObs module. 
    
`massRatioToMasses(msml, Nt, cbeta, Nf='21')`
 
    This is a way to get either the light+strange masses or quark+preconditioner masses
    given the ratio ms/ml or mpre/mf. This uses HotQCD ensembles, listed in tables below.

    Args:
        msml (int)
        Nt (int)
        cbeta (str): Beta as a character. e.g. 6.345 is 6345 
        Nf (str): Number of dynamical fermions. Defaults to '21'.

    Returns:
        str, str: either (ml, ms) or (mf, mpre) 
    
`quarkMassTableHISQ(Nf, Nt, msml) -> dict`

    Lookup tables for HotQCD HISQ runs.
    
`HotQCDParams(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf='21', scaleYear=None, mu=0)`

    A class to handle and check the input parameters of a lattice run, especially for HotQCD.
    
