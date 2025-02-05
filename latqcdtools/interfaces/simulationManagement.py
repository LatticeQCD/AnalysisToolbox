# 
# simulationManagement.py                                                               
# 
# D. Clarke
# 
# This is some Python code to help assist with managing lattice simulations. 
# 

import glob
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger

def countConfigurations(targetFolder,name,delimiter='.'):
    """
    Count the number of configurations in the target folder.

    Args:
        targetFolder (str)
        name (str): Assume configuration name has form name.###. 
        delimiter (str, optional): The delimiter between name and ###. Defaults to '.'.

    Returns:
        int: number of configurations in targetFolder 
    """
    checkType(str,targetFolder=targetFolder)
    checkType(str,name=name)
    checkType(str,delimiter=delimiter)
    confs=[]
    files=list(glob.iglob(f'{targetFolder}/{name}*'))
    for f in files:
        try:
            _, confno = f.rsplit(delimiter,1)
            confs.append(int(confno))
        except ValueError:
            logger.TBRaise(f'Configuration names in {targetFolder} must end in int. file={f}, delimiter="{delimiter}"')
    logger.info(f'{targetFolder}:  min={min(confs)}  max={max(confs)}')
    return len(confs)

