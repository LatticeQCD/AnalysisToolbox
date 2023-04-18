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


# --------------------- BASE TESTS


cd base
runTestRoutine testDataCleaner.py
rm Toolbox.log
runTestRoutine testLogger.py
_compareFiles control.log Toolbox.log
_checkPassError $? "testLogger" 
runTestRoutine testPrintErrorBars.py
runTestRoutine "testReadWrite.py --type str"
runTestRoutine "testReadWrite.py --type list"
runTestRoutine testUtilities.py
cd ..


# ---------------- INTERFACE TESTS


cd interfaces
runTestRoutine testReadWriteConf.py
cd ..


# --------------------- MATH TESTS


cd math
runTestRoutine testMath.py
runTestRoutine testDeriv.py
runTestRoutine testInt.py
runTestRoutine testPolynomial.py
runTestRoutine testSpline.py
runTestRoutine testSU3.py
cd ..


# --------------- STATISTICS TESTS


cd statistics
runTestRoutine testAutocor.py 
runTestRoutine testBootstrap.py 
runTestRoutine testErrorProp.py
runTestRoutine testFit.py
runTestRoutine testJackknife.py 
runTestRoutine testStats.py 
runTestRoutine testWeightedMean.py
cd ..


# ------------------ PHYSICS TESTS


cd physics
runTestRoutine testContExtrap.py
runTestRoutine testDenseObs.py
runTestRoutine testLatticeParams.py
runTestRoutine testPolyakovTools.py
runTestRoutine testScales.py
runTestRoutine testStaticPotential.py
runTestRoutine testStatPhys.py
runTestRoutine testHRG.py
cd ..


echo
echo "Done!"
echo
