#!/bin/bash

# 
# testingTools.bash
# 
# D. Clarke 
# 
# Some useful tools and definitions for carrying out LatticeToolbox tests. 
#

cgrn="\e[92m" 
cred="\e[91m"
cyan="\e[36m"
endc="\e[0m"

showOutputOnScreen=true

function runTestRoutine {
    routine="$1"
    outFile='OUT_'${routine%%.py*}
    errFile='runERR_'${routine%%.py*}
    echo
    echo '  '${routine}
    if [ ${showOutputOnScreen} ]; then
        echo "  -------------------------"
        echo
        python3 ./${routine}
    else
        python3 ./${routine} >> ${outFile} 2>> ${errFile}
        # This is just to remove these files if they are empty.
        if [ ! -s ${errFile} ]; then rm ${errFile}; fi
        if [ ! -s ${outFile} ]; then rm ${outFile}; fi
    fi
}


function _bashError {
  echo
  echo -e "  ${cred}ERROR: $1 ${endc}"
  echo
} 


function _bashPass {
  echo
  echo -e "  ${cgrn}SUCCESS: $1 ${endc}"
  echo
}


function _checkPassError {
  if [ $1 -eq 0 ]; then
    _bashPass "$2"
  else
    _bashError "$2"
  fi
}


function _compareFiles {
  diff $1 $2 &>/dev/null
  return $?
}
