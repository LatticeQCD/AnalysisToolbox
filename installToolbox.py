# 
# installToolbox.py
# 
# D. Clarke 
# 
# Add the AnalysisToolbox location to the PYTHONPATH in your bashrc.
#

import os

TOOLBOXDIR = os.getcwd()
BASHRCFILE = os.path.expanduser("~")+"/.bashrc"

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


