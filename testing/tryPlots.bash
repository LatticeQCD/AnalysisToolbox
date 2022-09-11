#!/bin/bash

# 
# tryPlots.bash                                                               
# 
# D. Clarke 
# 
# Script to try all the plotting tests to see how it looks. 
#

source "testingTools.bash"

echo
echo "Running plotting tests:"
echo

cd plotting 
runTestRoutine plot_random.py
runTestRoutine test_plot_formats.py
runTestRoutine test_plot.py

echo
echo "Done!"
echo
