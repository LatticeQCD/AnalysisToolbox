#!/bin/bash

# 
# runTests.bash                                                               
# 
# D. Clarke 
# 
# Script to run every test in the AnalysisToolBox. Please add your tests to this script
# when you have written them. 
#

source "testingTools.bash"

echo
echo "Running python tests:"
echo


# --------------- STATISTICS TESTS


cd statistics
runTestRoutine gaudifTest.py
runTestRoutine testautocor.py 
runTestRoutine testBootstrap.py 
runTestRoutine testJackknife.py 
cd ..


# ------------------ PHYSICS TESTS


cd physics
runTestRoutine testLatticeParams.py
runTestRoutine testPolyakovTools.py
runTestRoutine testScales.py
cd ..


# ------------------ FITTING TESTS

cd fitting
runTestRoutine test_fit.py
cd ..

echo
echo "Done!"
echo
