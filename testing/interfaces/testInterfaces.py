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

    ref = latexTable() 
    ref.append(['1','2','3'])
    ref.append(['1','2','3'])
    ref.append(['1','2','3'])
    ref.outputTable('testInterface.tex')
    test = latexTable()
    test.readTable('testInterface.tex')
    lpass *= test==ref

    ref = redmineTable() 
    ref.append(['1','2','3'])
    ref.append(['1','2','3'])
    ref.append(['1','2','3'])
    ref.outputTable('testInterface.redmine')
    test = redmineTable()
    test.readTable('testInterface.redmine')
    lpass *= test==ref

    writeYAML(refDict,'test.yaml')
    testDict = readYAML('test.yaml')
    lpass *= testDict==refDict

    deleteFile('test.yaml') 
    deleteFile('testInterface.tex') 
    deleteFile('testInterface.redmine') 

    concludeTest(lpass) 

if __name__ == '__main__':
    testInterfaces()