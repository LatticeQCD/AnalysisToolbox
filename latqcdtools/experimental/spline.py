import numpy as np
import math
from scipy.linalg import inv

"""
Evaluate nth derivative of (x - knot)**order

Parameters
----------
x:
    position to evaluate
knot:
    shift of the polynomial
order:
    Order of the polynomial
nderiv:
    nth derivative
"""

def pol_deriv(x, knot, order, nderiv):
    if nderiv > order:
        return 0
    return math.factorial(order) / math.factorial(order - nderiv) * (x - knot)**(order - nderiv)
    

def jump_func(x, knot, order, nderiv):
    return (x >= knot) * pol_deriv(x, knot, order, nderiv)


"""
Simple spline.
Parameters
----------

x: scalar or array_like
    Position where to evaluate the spline
knots: array_like
    Positions of the knots
coeffs: array_like
    coefficients of the spline. 
base_point: scalar
    Position at which base polynomial is evaluated
order: int
    order of the spline
nsteady_deriv: int
    number of derivatives that should be steady. Default = order - 1
nderiv: int
    order of the derivative to be computed. Default = 0
"""

def spline(x, coeffs, knots, order = 2, nsteady_deriv = None, nderiv = 0, base_point = 0):

    if nsteady_deriv is None:
        nsteady_deriv = order - 1

    nnot_steady = order - nsteady_deriv

    ncoeffs = nnot_steady * len(knots) + order + 1

    if len(coeffs) != ncoeffs:
        raise ValueError("Number of coefficients does not match to order and number of knots")

    ret = 0

    for i in range(order + 1):
        ret += coeffs[i] * pol_deriv(x, base_point, i, nderiv)

    for j in range(0, len(knots)):
        for k in range(nnot_steady):
            ret += coeffs[order + 1 + j * nnot_steady + k] * jump_func(x, knots[j],
                    order - (nnot_steady - k - 1), nderiv)

    return ret


"""Calculate the number of spline coefficients based on spline parameters"""

def get_num_spline_coeffs(order, nsteady_deriv, nknots, constraints):

    if nsteady_deriv is None:
        nsteady_deriv = order - 1

    nnot_steady = order - nsteady_deriv
    return order + 1 + nnot_steady * nknots - len(constraints)



"""
Constrain a linear function

Parameters
----------
func: Callable
    Function that shall be constraints. It has to look like
    func(x, params, *args, **kwargs, nderiv = 0)
    where x is the position at which the function can be evaluated. The parameters that enter
    into the function are passed with params. *args and **kwargs are optional parameters.
    nderiv is a flag to calculate the n'th derivative of the function. This is necessary for 
    constraints that work on derivatives.

params: array_like
    parameters that should not be constraints. This means that the parameters that are used to 
    fulfill the constraints should not be included in this array. They will be put in.

index_sel: Callable
    Function that calculates the index of the parameters that are used for the constraints.
    The indices are the indices of the full parameter array.

constraints:
    Constraints which the function has to follow:
    [[x-position of constraint, nth-derivative, value], [...]]

func_args:
    optional args for func

func_kwargs:
    optional kwargs for func

Returns
-------

Array of new coefficients for the function
"""

def constraint_linear_function(func, params, index_sel, constraints, func_args = (),
        func_kwargs = {}):

    params = np.array(params, dtype = float)

    nconstrains = len(constraints)

    param_inds = index_sel(params, constraints, *func_args, **func_kwargs)

    for ind in sorted(param_inds):
        params = np.insert(params, ind, 0)

    
    rhs = np.empty(nconstrains)
    mat = np.empty(shape = (nconstrains, nconstrains))

    for i, constraint in enumerate(constraints):

        pos = constraint[0]
        nderiv = constraint[1]
        val = constraint[2]

        for j, (ind, constraint) in enumerate(zip(param_inds, constraints)):

            contrib_params = np.zeros_like(params)
            contrib_params[ind] = 1

            mat[i][j] = func(pos, contrib_params, *func_args, nderiv = nderiv, **func_kwargs)

        cur_func_value = func(pos, params, *func_args, nderiv = nderiv, **func_kwargs)
        rhs[i] = val - cur_func_value

    if len(constraints) > 0:
        constraint_params = inv(mat).dot(rhs)

        for ind, param in zip(param_inds, constraint_params):
            params[ind] = param

    return params


"""Index selection function for simple splines"""

def spline_ind_sel(params, constraints, knots, order, nsteady_deriv, **kwargs):

    nnot_steady = order - nsteady_deriv

    indlist = []

    for constraint in constraints:

        relevant_knot_pos = np.searchsorted(knots, constraint[0]) - 1
        ind = (relevant_knot_pos + 1) * nnot_steady + order

        ind_in_list = True

        while ind_in_list:
            ind_in_list = ind in indlist

            if ind_in_list:
                ind -= 1
                if ind < 0:
                    raise IndexError("No parameter left to fulfill "
                            "constraint + str(constraint)")
            else:
                indlist.append(ind)

    return indlist



"""
Constraint simple spline. For each constraint, one coefficient is computed based on that
constraint. The syntax for the constraints is as follows:
    [[x-position of constraint, nth-derivative, value], [...]]
For each knot interval, it is checked whether a constraint determines the spline
coefficient in that knot interval. If so, the coefficient is calculated based on the
constrained value.

Examples
--------
>>> constraint_spline(1.5, [1,4], [1,2,5], order = 1, constraints=[[1.5, 0, 1]])
1.0

The missing coefficient is placed at the corresponding knot position, so that the final
coefficients evaluate to [1, 2, -6.0, 5]

If there is more than one coefficient available, the latter ones are used.
This is especially the case if there is constraints before the first knot, or if nsteady_deriv
is smaller than order - 1:

>>> constraint_spline(0.5, [1,4], [2,2,5], order = 2, constraints=[[0.5, 0, 1], [0.5, 1, 0]])
1.0
>>> constraint_spline(0.5, [1,4], [2,2,5], order = 2, constraints=[[0.5, 0, 1], [0.5, 1, 0]],
>>>        nderiv = 1)
0.0

The used coefficients are [2, -4.0, 4.0, 2, 5]



Parameters
----------

x: scalar or array_like
    Position where to evaluate the spline
knots: array_like
    Positions of the knots
coeffs: array_like
    coefficients of the spline without the coeffecients that are computed from the contraints
base_point: scalar
    Position at which base polynomial is evaluated
order: int
    order of the spline
nsteady_deriv: int
    number of derivatives that should be steady. Default = order - 1
nderiv: int
    order of the derivative to be computed. Default = 0

constraints: 2D array_like
    See above text
"""

def constraint_spline(x, knots, coeffs, order = 2, nsteady_deriv = None,
        nderiv = 0, constraints = [], base_point = 0):

    if nsteady_deriv is None:
        nsteady_deriv = order - 1

    for constraint in constraints:
        if constraint[0] == base_point:
            raise ValueError("Constraint at base_point = "
                    + str(base_point) + ". Please shift base_point")


    if len(coeffs) != get_num_spline_coeffs(order, nsteady_deriv, len(knots), constraints):
        raise ValueError("Number of coefficients does not match to number"
                " of knots, order and constraints")


    new_coeffs = constraint_linear_function(spline, coeffs, spline_ind_sel, constraints,
            func_args = (knots,),
            func_kwargs = {'order' : order,
                'nsteady_deriv' : nsteady_deriv,
                'base_point': base_point
                } )


    return spline(x, new_coeffs, knots, order, nsteady_deriv, nderiv, base_point = base_point)
