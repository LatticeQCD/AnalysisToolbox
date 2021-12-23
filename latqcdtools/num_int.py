from numpy.polynomial.laguerre import laggauss
import numpy as np

deg_cache = 0
x_cache = None
w_cahce = None

def int_lag_gauss(func, args = (), deg = 100):
    x, w = laggauss(deg)
    res = 0
    for i in range(deg):
        res += func(x[i], *args) * w[i]
    return res


def epsilon(data):
    N = len(data)

    # Check if we get an array of values
    try:
        data[0][0]
        e = np.zeros((N,N,data.shape[1]))
    except IndexError:
        e = np.zeros((N,N))

    for i in range(N):
        e[i][0] = data[i]
    
    for k in range(0, N - 1):
        e[k][1] = 1 / (e[k+1][0] - e[k][0])
    

    for p in range(2, N):
        for k in range(0, N - p):
            e[k][p] = e[k+1][p-2] + 1 / (e[k+1][p-1] - e[k][p-1])
        
    if (N - 1) % 2 == 0:
        return e[0][N-1]
    else:
        return e[1][N-2]