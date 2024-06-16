# 
# autoGlossary.py                                                               
# 
# D. Clarke 
# 
# Go through latqcdtools and generate a .md file for each module. That .md file
# should contain all classes and functions along with their docstrings. All
# .md files are saved under docs_src/glossary to be incorporated into the
# documentation easily. 
#

import os, pkgutil, importlib, inspect

LIBRARY_NAME = 'latqcdtools'
LIBRARY_PATH = os.path.dirname(importlib.import_module(LIBRARY_NAME).__file__)


def get_functions_with_docs(module):
    functions = inspect.getmembers(module, inspect.isfunction)
    module_name = module.__name__
    functions = [(name, func) for name, func in functions if func.__module__ == module_name]
    return [(name, inspect.signature(func), func.__doc__) for name, func in functions]

def get_classes_with_docs(module):
    functions = inspect.getmembers(module, inspect.isclass)
    module_name = module.__name__
    functions = [(name, func) for name, func in functions if func.__module__ == module_name]
    try:
        return [(name, inspect.signature(func), func.__doc__) for name, func in functions]
    except ValueError:
        return [(name, None, func.__doc__) for name, func in functions]

glossary = open('docs_src/glossary/glossary.md','w')

glossary.write("Glossary\n")
glossary.write("=============\n")
glossary.write("\n")
glossary.write("Here is a glossary of all classes and methods in the AnalysisToolbox organized by module.\n")
glossary.write("\n")
glossary.write("```{toctree}\n")
glossary.write("---\n")
glossary.write("maxdepth: 1\n")
glossary.write("---\n")

for importer, module_name, ispkg in pkgutil.walk_packages(path=[LIBRARY_PATH], prefix=LIBRARY_NAME + '.'):

    if not ispkg:
        module = importlib.import_module(module_name)

        glossary.write(module_name+".md\n")

        mdfile = open('docs_src/glossary/'+module_name+'.md','w') 

        mdfile.write(module_name+'\n')
        mdfile.write('=============\n\n')

        for func_name, sig, func_doc in get_functions_with_docs(module):
            doc = func_doc
            if func_doc is None:
                doc = "\n"
            mdfile.write(f"`{func_name}{sig}`\n{doc}\n")

        for func_name, func_sig, func_doc in get_classes_with_docs(module):
            doc = func_doc
            if func_doc is None:
                doc = "\n"
            sig = func_sig
            if func_sig is None:
                sig = "\n"
            mdfile.write(f"`{func_name}{sig}`\n{doc}\n")
        
        mdfile.close()

glossary.write("```\n")
glossary.close()