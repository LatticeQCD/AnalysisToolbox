#!/bin/bash

# 
# testingTools.bash
# 
# D. Clarke 
# 
# Some useful tools and definitions for carrying out AnalysisToolBox tests. 
#

cred="\e[91m"
cyan="\e[36m"
endc="\e[0m"

showOutputOnScreen=true

function runTestRoutine {
    routine="$1"
    outFile='OUT_'${routine%%.py*}
    errFile='runERR_'${routine%%.py*}
    echo '  '${routine}
    if [ ${showOutputOnScreen} ]; then
        echo "  -------------------------"
        python ./${routine}
        echo
    else
        python ./${routine} >> ${outFile} 2>> ${errFile}
        # This is just to remove these files if they are empty.
        if [ ! -s ${errFile} ]; then rm ${errFile}; fi
        if [ ! -s ${outFile} ]; then rm ${outFile}; fi
    fi
}
