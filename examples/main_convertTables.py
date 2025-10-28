# 
# main_convertTables.py                                                               
# 
# D. Clarke
# 
# An example showing how the AnalysisToolbox can help with reading LaTeX tables
# and converting them to something usable by python. Note the AnalysisToolbox
# also has some helpful structures for csv and Wikipedia tables. 
# 
import re
from latqcdtools.interfaces.interfaces import latexTable
import latqcdtools.base.logger as logger

# This is an example table from Phys. Lett. B 858 (2024) 139040
table1 = latexTable() 
table1.readTable('../datasets/example_files/table1.tex')

#
# What follows are some functions to convert from LaTeX symbols to
# math strings. These are not coded in the AnalysisToolbox because
# LaTeX tables can vary wildly, which means it's better for the user
# to decide how to parse.
#
def removeBF(latex_string):
    match = re.search(r'\$\\mathbf\{(.+?)\}\$', latex_string)
    if match:
        return match.group(1)
    return latex_string 

def convert_minus_fraction(latex_string):
    pattern = r'-\\frac\{(.+?)\}\{(.+?)\}'
    match = re.search(pattern, latex_string)
    if match:
        numerator = match.group(1)
        denominator = match.group(2)
        return f"-{numerator}/{denominator}"
    return latex_string

def convert_fraction(latex_string):
    pattern = r'\\frac\{(.+?)\}\{(.+?)\}'
    match = re.search(pattern, latex_string)
    if match:
        numerator = match.group(1)
        denominator = match.group(2)
        return f"{numerator}/{denominator}"
    return latex_string

def removeMM(latex_string):
    pattern = r'\$(.+?)\$'
    match = re.search(pattern, latex_string)
    if match:
        return match.group(1)
    return latex_string

#
# Using the above functions, let's clean up our table. We can see now that
# each list newRow is much more easily used by Python. 
#
for row in table1:
    newRow = []
    for ele in row:
        newRow.append(removeMM(removeBF(convert_fraction(convert_minus_fraction(ele)))))
    logger.info(newRow)
    
