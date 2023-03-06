# 
# testInt.py                                                               
# 
# D. Clarke
# 
# Quick test of the integration class.
# 


from latqcdtools.math.num_int import integrateData, integrateFunction
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def f(s):
    return s


def testInt():

    x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
         30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48]
    y = [ 1.67787710e+01,  1.25847027e+01,  1.10282977e+01,  5.40264891e+00,
          3.12776222e+00, -4.17877810e-01, -1.56451625e+00, -3.05995384e+00,
         -3.25481867e+00, -3.58460190e+00, -3.33892230e+00, -3.18112342e+00,
         -2.77937463e+00, -2.49030668e+00, -2.10789478e+00, -1.80839245e+00,
         -1.49496561e+00, -1.27000897e+00, -1.05991815e+00, -8.91894201e-01,
         -6.92139494e-01, -5.52904745e-01, -4.16482841e-01, -3.46622975e-01,
         -2.26077138e-01, -2.13762415e-01, -1.13309514e-01, -1.47992939e-01,
         -7.79446646e-02, -1.26137554e-01, -7.30776714e-02, -7.82699808e-02,
         -1.52375973e-02, -5.90653747e-03,  6.46070820e-02,  8.73126699e-02,
          1.09933589e-01,  1.15683939e-01,  7.51285569e-02,  7.17124076e-02,
          3.87964264e-02,  3.71900954e-02,  3.40575679e-02,  3.04882265e-02,
          2.01035478e-03, -2.99829433e-02, -2.63919195e-02, -3.27346097e-02,
         -4.18786374e-02]

    EPSILON =1e-8
    I_simp  =3.86033230
    I_trap  =5.69923345
    I_spline=5.49427272

    # a result from real lattice data
    print_results(integrateData(x,y,method='simpson')  ,I_simp  ,text="simpson"  ,prec=EPSILON)
    print_results(integrateData(x,y,method='trapezoid'),I_trap  ,text="trapezoid",prec=EPSILON)
    print_results(integrateData(x,y,method='spline')   ,I_spline,text="spline"   ,prec=EPSILON)

    # try a simple integral to compare against the fundamental theorem of calculus
    FTC=(x[-1]**2-x[0]**2)/2

    print_results(integrateData(x,x,method='simpson')  ,FTC,text="simpson integrate f(x) = x")
    print_results(integrateData(x,x,method='trapezoid'),FTC,text="trapezoid integrate f(x) = x")
    print_results(integrateData(x,x,method='spline')   ,FTC,text="spline integrate f(x) = x")

    # try some vectorization
    print_results(integrateFunction(f,[0,0],[1,2],method='trapezoid'), [0.5,2], text="trapezoid vector function")
    print_results(integrateFunction(f,0,1,method='trapezoid'), 0.5, text="trapezoid scalar function")
    print_results(integrateFunction(f,[0,0],[1,2],method='quad'), [0.5,2], text="quadrature vector function")
    print_results(integrateFunction(f,0,1,method='quad'), 0.5, text="quadrature scalar function")


if __name__ == '__main__':
    testInt()