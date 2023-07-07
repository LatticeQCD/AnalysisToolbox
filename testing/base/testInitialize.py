# 
# testInitialize.py                                                               
# 
# D. Clarke 
# 
# Test the initializer. 
#

from latqcdtools.base.initialize import initialize, finalize

def testInitialize():
    initialize()
    finalize()

if __name__ == '__main__':
    testInitialize()
