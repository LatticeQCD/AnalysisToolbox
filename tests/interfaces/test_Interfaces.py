# 
# testInterfaces.py                                                               
# 
# D. Clarke
# 
# Testing of some generic interfacing tools. 
# 

from latqcdtools.interfaces.interfaces import latexTable, redmineTable, readYAML, writeYAML, \
    writeJSON, readJSON, convertTable, csvTable, paramFromEnsLabel
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
    ref.writeTable('testInterface.tex')
    test = latexTable()
    test.readTable('testInterface.tex')
    lpass *= test==ref

    ref = redmineTable() 
    ref.append(['1','2','3'])
    ref.append(['1','2','3'])
    ref.append(['1','2','3'])
    ref.writeTable('testInterface.redmine')
    test = redmineTable()
    test.readTable('testInterface.redmine')
    lpass *= test==ref

    convertTable('testInterface.tex','testInterface.csv',targetDelimiter=',')
    test = csvTable(delimiter=',')
    test.readTable('testInterface.csv')
    lpass *= test==ref

    writeYAML(refDict,'test.yaml')
    testDict = readYAML('test.yaml')
    lpass *= testDict==refDict

    writeJSON(refDict,'test.json')
    testDict = readJSON('test.json')
    lpass *= testDict==refDict

    deleteFile('test.yaml') 
    deleteFile('test.json') 
    deleteFile('testInterface.tex') 
    deleteFile('testInterface.redmine') 
    deleteFile('testInterface.csv')

    param = paramFromEnsLabel('l248f111b37000m00139736m00291117m0585145',format='MILC')
    if param != (24, 8, '111', '37000', '00139736', '00291117', '0585145'):
        logger.TBFail('paramFromEnsLabel')
        lpass = False

    concludeTest(lpass) 

if __name__ == '__main__':
    testInterfaces()