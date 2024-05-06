#!/bin/bash

# 
# runTests.bash                                                               
# 
# D. Clarke 
# 
# Script to run every test in the AnalysisToolbox. Please add your tests to this script
# when you have written them. 
#

source "testingTools.bash"

echo
echo "Running python tests:"
echo


# --------------------- BASE TESTS


cd base
runTestRoutine testCheck.py
runTestRoutine testDataCleaner.py
runTestRoutine testInitialize.py
runTestRoutine testLogger.py
runTestRoutine testPrintErrorBars.py
runTestRoutine testReadWrite.py
runTestRoutine testSpeedify.py
runTestRoutine testUtilities.py
cd ..


# ---------------- INTERFACE TESTS


cd interfaces
runTestRoutine testReadWriteConf.py
runTestRoutine testInterfaces.py
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
runTestRoutine testBMA.py 
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
#runTestRoutine testDenseObs.py
runTestRoutine testGauge.py
runTestRoutine testHotQCDEos.py
runTestRoutine testHRG.py
runTestRoutine testIdeal.py
runTestRoutine testLatticeParams.py
runTestRoutine testPolyakovTools.py
runTestRoutine testRunningCoupling.py
runTestRoutine testScales.py
runTestRoutine testStaticPotential.py
runTestRoutine testStatPhys.py
runTestRoutine testUnits.py
cd ..


runTestRoutine testLegacy.py


echo
echo "Done!"
echo
