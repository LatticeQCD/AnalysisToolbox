
"""
Calculating the Taylor co-efficients for HRG. currently it calculates only pressure co-efficients
Create by Jishnu Goswami 19/07/2022
"""


import math as mp
import numpy as np
import sympy as sy
from sympy import Symbol,symbols,Indexed
from latqcdtools.physics.HRG import HRG,EV_HRG
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import getArgs
import re
#from sympy import init_printing
#init_printing()

#P=Symbol('P')

def string_search(arr):
    temp = [] 
    for t in arr:
        temp.append(re.findall("[0-9][0-9][0-9]", t)[0])
    return temp


muB,muQ,muS,s1,s3,s5,s7,s9,q1,q3,q5,q7,q9=symbols('muB,muQ,muS,s1,s3,s5,s7,s9,q1,q3,q5,q7,q9')


order=6
#P=sum(1/sy.factorial(i)*1/sy.factorial(j)*1/sy.factorial(k) * muB ** i * muQ ** j* muS ** k* symbols('X.{}.{}.{}'.format(i,j,k)) for i in range(order) for j in range(order) for k in range(order) if (i+j+k)%2==0 and i+j+k < order)
P=sum(1/sy.factorial(i)*1/sy.factorial(j)*1/sy.factorial(k) * muB ** i * muQ ** j* muS ** k* symbols('X{}{}{}'.format(i,j,k)) for i in range(order) for j in range(order) for k in range(order) if (i+j+k)%2==0 and i+j+k < order)

"""
Co-efficients of Pressure
"""
print ("Calculations of Co-efficients of Pressure....")
m = 4
f = P.subs({muS:s1*muB + s3*muB**3 + s5*muB**5 + s7*muB**7 + s9*muB**9 , muQ:q1*muB + q3*muB**3 + q5*muB**5 + q7*muB**7 + q9*muB**9})
F = f.expand()
PBm = mp.factorial(m)*F.coeff(muB**m)  # Coefficients will come without the factorial


test = str(PBm).split('+')

BQS = string_search(test)


hadrons,M,Q,B,S,C,g = np.loadtxt("../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))


obs = np.genfromtxt("s_q_QM_hadron_list_ext_strange_2020.txt",names = True)

T  = obs['TMeV']
s1 = obs['s1']
s3 = obs['s3']
q1 = obs['q1']
q3 = obs['q3']

#T = np.linspace(130, 180, 101)
w  = np.array([1 if ba==0 else -1 for ba in B])
QMhrg      = HRG(M,g,w,B,S,Q)

for i in range(len(test)):
    Border = int(BQS[i][0])
    Qorder = int(BQS[i][1])
    Sorder = int(BQS[i][2])
    globals()['X%s'%BQS[i]] = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder)

P2n = eval(str(PBm))

#print (P2n)

np.savetxt("pressure_%dorder_coeff_nS_zero.txt"%m,np.c_[T , P2n], fmt = '%0.6e')

