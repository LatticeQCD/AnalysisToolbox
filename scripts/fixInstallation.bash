# 
# fixInstallation.bash
# 
# D. Clarke 
# 
# It can happen that sometimes a local upgrade to Python3 will break
# some of the packages in requirements.txt; for example an upgrade can
# lead to a conflict between numba and numpy. This script can sometimes
# fix these problems. 
#

pip3 install --upgrade --force-reinstall -r ../requirements.txt
