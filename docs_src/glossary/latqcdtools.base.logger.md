latqcdtools.base.logger
=============

`TBError(*args, frame=2)`

    Print error message and exit with -1.

    Args:
        frame (int, optional): Controls the name of the caller. Defaults to method that called TBError.
    
`TBFail(*args)`


`TBPass(*args)`


`TBRaise(*args, frame=2)`

    Print error message and raise ToolboxException. 

    Args:
        frame (int, optional): Controls the name of the caller. Defaults to method that called TBRaise.
    
`_getCallerName(frame)`
 
    Gets the name of the function that calls the present function. 
    
`_getTimeStamp()`
 
    Get HH:MM:SS 
    
`_log(outString)`


`createLogFile(filename='Toolbox.log')`
 
    Have output sent also to a log file filename. If this file already exists, it will get deleted. We use the
    logging module because it knows how to handle multiple processes writing to the same file. 
    
`debug(*args, frame=2)`


`details(*args)`


`info(*args)`


`progress(*args)`


`set_log_level(level)`


`warn(*args, frame=2)`


`ToolboxException
`


