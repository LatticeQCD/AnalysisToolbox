#!/bin/python3

# 
# configureToolbox.py
# 
# D. Clarke 
# 
# Set up your AnalysisToolbox. Depending on your preferences, this either
#
#   1. Adds a line to your .bashrc
#   2. Sets up a venv for you with the name AnalysisToolbox 
#


# ----------------------------------------------------------------------- USER PREFERENCES

SHELL    = "BASH"    # POWERSHELL to be added later
STRATEGY = "BASIC"   # Options are "VENV" and "BASIC"

# ------------------------------------------------------------------------------ MAIN CODE 

import os
from subprocess import run, PIPE

def shell(*args):
    """ Carry out the passed arguments args in the shell. Can be passed as a single
        string or as a list. Captures and returns output of shell command. E.g.
            shell('ls -lah')
    """
    args = [str(s) for s in args]
    process = run(' '.join(args),shell=True,check=True,stdout=PIPE,universal_newlines=True)
    return process.stdout

TOOLBOXDIR = os.getcwd()
BASHRCFILE = os.path.expanduser("~")+"/.bashrc"

ALLOWED_STRATEGIES = ["VENV","BASIC"]
ALLOWED_SHELLS = ["BASH"]

if not SHELL in ALLOWED_SHELLS:
    print("Unrecognized shell",SHELL)
    exit(-1)

if not STRATEGY in ALLOWED_STRATEGIES:
    print("Unrecognized strategy",STRATEGY)
    exit(-1)

if STRATEGY=="VENV":

    shell('mkdir -p venv; cd venv; python3 -m venv AnalysisToolbox')
    

found = False

# Read the file and check if it contains a line starting with "export PYTHONPATH="
with open(BASHRCFILE, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if line.startswith("export PYTHONPATH="):
            found = True
            print("Appending AnalysisToolbox to PYTHONPATH.")
            lines[i] = line.strip()[:-1] + ':'+TOOLBOXDIR+'"\n'
            break

# If the line is not found, create a new line 
if not found:
    print("Adding PYTHONPATH line for the AnalysisToolbox.")
    lines.append('export PYTHONPATH="${PYTHONPATH}:'+TOOLBOXDIR+'"\n')

# Write the updated content back to the file
with open(BASHRCFILE, 'w') as file:
    file.writelines(lines)
