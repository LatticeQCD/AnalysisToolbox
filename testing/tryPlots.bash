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

cd base 
runTestRoutine plotRandom.py
runTestRoutine testPlotFormats.py
runTestRoutine testPlot.py
runTestRoutine testHist.py

echo
echo "Done!"
echo
