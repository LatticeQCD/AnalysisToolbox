# 
# testInterfaces.py                                                               
# 
# D. Clarke
# 
# Testing of some generic interfacing tools. 
# 

from latqcdtools.interfaces.interfaces import latexTable, redmineTable, readYAML, writeYAML 
from latqcdtools.testing import concludeTest
from latqcdtools.base.utilities import deleteFile
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

refDict = { 'amuL':
    { 'priors' :
        { '[0, 0.4]' :
            { 'a2' :
                { 'Raw' :
                    { 'Ca2': '0(10000)',
                      'Cs': '0.0(3)'
                    }
                },
              'a4':
                { 'Raw':
                    { 'Ca2': '0(10000)',
                      'Ca4': '0(10000)'
                    }
                },
              'a6':
                { 'Raw':
                    { 'Ca2': '0(10000)',
                      'Ca4': '0(10000)'
                    }
                }
            }
        }
    }
}


def testInterfaces():

    lpass = True

    test = latexTable() 
    test.append([1,2,3])
    test.append([1,2,3])
    test.append([1,2,3])
    test.outputTable('testInterface.tex')

    test = redmineTable() 
    test.append([1,2,3])
    test.append([1,2,3])
    test.append([1,2,3])
    test.outputTable('testInterface.redmine')

    # Test identity operation
    writeYAML(refDict,'test.yaml')
    testDict = readYAML('test.yaml')
    lpass *= testDict==refDict
    deleteFile('test.yaml') 

    concludeTest(lpass) 

if __name__ == '__main__':
    testInterfaces()