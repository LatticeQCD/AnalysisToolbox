latqcdtools.base.readWrite
=============

`_lab(num) -> str`
 
    Create a short string label for each column of a data table. Needed for writeTable. 
    
`readTable(filename, unpack=True, col=None, minVal=-inf, maxVal=inf, **kwargs) -> numpy.ndarray`
 
    Wrapper for np.loadtxt. It unpacks by default to prevent transposition errors, and also optionally
    allows the user to restrict the table based on the range of one of the columns.

    Args:
        filename (str):
            Name of the file to open. 
        unpack (bool, optional): 
            If False, reads the table in a confusing, useless way. Defaults to True.
        col (int, optional): 
            Use this column to restrict the range of the rest of the table. Defaults to None.
        minVal (float, optional): 
            Minimum value for above restriction. Defaults to None.
        maxVal (float, optional):
            Maximum value for above restriction. Defaults to None.

    Returns:
        np.array: Data table. 
    
`writeTable(filename, *args, **kwargs)`
 
    Wrapper for np.savetxt. The idea is that you can use it like this:
    
    writeTable('file.txt',col1,col2,header=['header1','header2])
    
    This works for an arbitrary number of 1-d columns col. It seems much more intuitive to me 
    that you pass columns as arguments than whatever np.savetxt is doing.

    Args:
        filename (str): output file name
    
