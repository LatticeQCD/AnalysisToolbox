# 
# testReadWrite.py                                                               
# 
# D. Clarke
# 
# Tests for read and write methods.
#

import numpy as np
from latqcdtools.base.readWrite import readTable, writeTable
from latqcdtools.testing import print_results, concludeTest


def testReadWrite():

    lpass = True

    xdata = np.linspace(0,1,101)
    ydata = xdata**2
    zdata = ydata + np.pi*1j

    longString="Mitten_in_dem_kleinen_Teiche_steht_ein_Pavillon_aus_grünem_und_aus_weißem_Porzellan"

    textColumn=[longString for i in range(len(xdata))]

    writeTable("test.d",textColumn,xdata,ydata,zdata,header=["long string","xdata","ydata","zdata"],width=20,digits=10)
    stest,xtest,ytest,Reztest,Imztest = readTable("test.d",dtype='U99,f8,f8,f8,f8')

    ztest = Reztest + 1j*Imztest

    lpass *= print_results(xdata,xtest,text='x data')
    lpass *= print_results(ydata,ytest,text='y data')
    lpass *= print_results(zdata,ztest,text='z data')

    concludeTest(lpass)

if __name__ == '__main__':
    testReadWrite()