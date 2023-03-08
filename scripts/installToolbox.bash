#!/bin/bash

cd ..
currentDirectory=$(pwd)
bashrcFile=${HOME}/.bashrc

echo "" >> ${bashrcFile}
echo "# For the AnalysisToolbox" >> ${bashrcFile}
echo "export PYTHONPATH=\"\${PYTHONPATH}:${currentDirectory}\"" >> ${bashrcFile}

source ${bashrcFile}

echo
echo "Done! If you are using the AnalysisToolbox locally," 
echo "you can now simply"
echo "  pip3 install -r requirements.txt"
echo "Otherwise if you are on a supercomputer, you may need"
echo "to check its documentation to see which modules must"
echo "be loaded."
echo
