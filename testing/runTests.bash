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
runTestRoutine "testReadWrite.py --type str"
diff test.d control.d
runTestRoutine "testReadWrite.py --type list"
diff test.d control.d
runTestRoutine testUtilities.py
runTestRoutine testPrintErrorBars.py
cd ..


# --------------------- MATH TESTS


cd math
runTestRoutine testMath.py
runTestRoutine testDeriv.py
runTestRoutine testInt.py
runTestRoutine testPolynomial.py
runTestRoutine testSpline.py
cd ..


# --------------- STATISTICS TESTS


cd statistics
runTestRoutine testGauDif.py
runTestRoutine testautocor.py 
runTestRoutine testBootstrap.py 
runTestRoutine testJackknife.py 
runTestRoutine testErrorProp.py
runTestRoutine testFit.py
runTestRoutine testWeightedMean.py
cd ..


# ------------------ PHYSICS TESTS


cd physics
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
