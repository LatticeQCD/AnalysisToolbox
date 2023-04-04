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
runTestRoutine plotRandom.py
runTestRoutine testPlotFormats.py
runTestRoutine testPlot.py

echo
echo "Done!"
echo
