def SIGN(a, b):
    if b >= 0.0:
        return abs(a)
    else:
        return abs(-a)


def solve(func, y, x1, x2, tol, *args):
    def func_zero(x):
        return func(x, *args) - y

    ITMAX = 1000
    EPSILON = 3.0e-14
    a = x1
    b = x2
    c = x2
    fa = func_zero(a)
    fb = func_zero(b)
    if (fa > 0.0 and fb > 0.0) or (fa < 0.0 and fb < 0.0):
        raise ValueError("Root must be bracketed in zbrenz")

    fc = fb
    for itt in range(0, ITMAX):
        if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
            c = a
            fc = fa
            e = d = b - a
        if abs(fc) < abs(fb):
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        tol1 = 2.0 * EPSILON * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        if abs(xm) <= tol1 or fb == 0.0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            min1 = 3.0 * xm * q - abs(tol1 * q)
            min2 = abs(e * q)
            if 2.0 * p < min1 and 2.0 * p < min2:
                e = d
                d = p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d
        a = b
        fa = fb
        if abs(d) > tol1:
            b += d
        else:
            b += SIGN(tol1, xm)
        fb = func_zero(b)
    print("Maximum number of iterations exceeded in zbrent")
    return 0