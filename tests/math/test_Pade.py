#
# test_Pade.py
#
# D. Clarke
#
# Testing for the singlePointPade class.
#

import numpy as np
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.math.pade import singlePointPade


def _series(num, den, N):
    """
    Taylor coefficients c[0..N-1] about x=0 of the rational function
    (sum num[i] x**i) / (sum den[j] x**j), assuming den[0] != 0.
    """
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    c = np.zeros(N)
    for k in range(N):
        acc = num[k] if k < len(num) else 0.0
        for j in range(1, min(k, len(den) - 1) + 1):
            acc -= den[j] * c[k - j]
        c[k] = acc / den[0]
    return c


def testPade():

    lpass = True

    # Taylor coefficients of exp about 0.
    c_exp = np.array([1.0, 1.0, 1/2, 1/6, 1/24, 1/120, 1/720])

    # [1/1] Padé of exp is (1 + x/2) / (1 - x/2).
    R11 = singlePointPade(c_exp, p=1, q=1)
    lpass *= print_results(R11.num_coeffs, [1.0,  0.5], text="exp [1/1] numerator",   prec=1e-12)
    lpass *= print_results(R11.den_coeffs, [1.0, -0.5], text="exp [1/1] denominator", prec=1e-12)

    # [3/3] Padé of exp should be very accurate near the origin.
    R33 = singlePointPade(c_exp, p=3, q=3)
    lpass *= print_results(R33(0.5), np.exp(0.5), text="exp [3/3] at x=0.5", prec=1e-6)

    # A genuinely rational function must be recovered exactly.
    num, den = [1.0, 2.0], [1.0, 3.0, 1.0]      # (1 + 2x) / (1 + 3x + x^2)
    c_rat = _series(num, den, 8)
    Rrat = singlePointPade(c_rat, p=1, q=2)
    xs = np.array([-0.3, 0.1, 0.7, 2.0])
    exact = (num[0] + num[1]*xs) / (den[0] + den[1]*xs + den[2]*xs**2)
    lpass *= print_results(Rrat(xs), exact, text="recover exact rational", prec=1e-10)

    # Expansion about x0 != 0. f(x) = 1/(1-x) about x0 = 2: coeffs are -(-1)**k.
    c_shift = np.array([-(-1.0)**k for k in range(6)])
    Rshift = singlePointPade(c_shift, p=1, q=1, x0=2.0)
    xs = np.array([1.5, 2.5, 3.0])
    lpass *= print_results(Rshift(xs), 1/(1-xs), text="Padé about x0=2", prec=1e-10)

    # Scalar in, scalar out.
    lpass *= print_results(float(np.isscalar(R33(0.5))), 1.0, text="scalar input returns scalar", prec=1e-12)

    # p + q must not exceed the order of the supplied Taylor series.
    threw = False
    try:
        singlePointPade(c_exp, p=4, q=4)     # needs 9 coefficients, only 7 given
    except Exception:
        threw = True
    lpass *= print_results(float(threw), 1.0, text="rejects p+q above series order", prec=1e-12)

    concludeTest(lpass)


if __name__ == '__main__':
    testPade()
