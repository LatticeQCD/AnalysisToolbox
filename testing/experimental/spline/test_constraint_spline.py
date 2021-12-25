#!/usr/bin/env python3

from latqcdtools import spline
from latqcdtools.tools import print_results


def deriv(x, c, t, d, constraints, nsteady, deriv):
    return diff_deriv(x, spline, args = (c, t, d, constraints, nsteady, deriv))


order = 3
knots = [1, 4]
coeffs = [1, 2, 4]
nsteady = 1

constraints = [
        [-2 , 0, 123],
        [1 , 1, -1.213],
        [3 , 1, -1],
        [3.5 , 2, 0.123],
        [5 , 0, 7.123]
        ]

print_results(spline.spline_ind_sel(coeffs, constraints, knots, order, nsteady), [3, 2, 5, 4, 7],
        text = "Index selection test")



for constraint in constraints:
    print_results(spline.constraint_spline(constraint[0], knots, coeffs, order, nsteady,
            constraint[1], constraints), constraint[2],
            text = str(constraint[1]) + "th derivative of Spline(" + str(constraint[0]) + ")",
            prec = 1e-10)



print("Second set of parameters:")


order = 3
knots = [1, 4]
coeffs = [1, 2]
nsteady = 2

constraints = [
        [-3 , 0, 123],
        [-2.5 , 1, -1.213],
        [3.5 , 2, 0.123],
        [5 , 0, 7.123]
        ]



print_results(spline.spline_ind_sel(coeffs, constraints, knots, order, nsteady), [3, 2, 4, 5],
        text = "Index selection test")

for constraint in constraints:
    print_results(spline.constraint_spline(constraint[0], knots, coeffs, order, nsteady,
            constraint[1], constraints), constraint[2],
            text = str(constraint[1]) + "th derivative of Spline(" + str(constraint[0]) + ")",
            prec = 1e-10)


print("Third set of parameters:")

order = 3
knots = [1, 4]
coeffs = [1, 2, 4]
nsteady = 1

constraints = [
        [0 , 0, 123],
        [1 , 1, -1.213],
        [4.1 , 1, -1],
        [4.5 , 2, 0.123],
        [5 , 0, 7.123]
        ]

print_results(spline.spline_ind_sel(coeffs, constraints, knots, order, nsteady), [3, 2, 7, 6, 5],
        text = "Index selection test")



for constraint in constraints:
    print_results(spline.constraint_spline(constraint[0], knots, coeffs, order, nsteady,
            constraint[1], constraints, base_point = -1), constraint[2],
            text = str(constraint[1]) + "th derivative of Spline(" + str(constraint[0]) + ")",
            prec = 1e-10)

