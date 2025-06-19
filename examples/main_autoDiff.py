# 
# main_autoDiff.py                                                               
# 
# D. Clarke
# 
# This is an instructional example where automatic differentiation is implemented. We
# compare the output against the derivative methods in the AnalysisToolbox. The example
# follows the pseudocode of https://en.wikipedia.org/wiki/Automatic_differentiation
#  

from latqcdtools.math.num_deriv import diff_deriv
import latqcdtools.base.logger as logger

# We need our objects to have .value and .partial methods to keep
# track of the value and partial derivative as we go.
class ValueAndPartial:
    def __init__(self, value, partial):
        self.value   = value
        self.partial = partial

class Expression:
    # We are going to automatically differentiate some general
    # expression. Each expression should have this method.
    def evaluate_and_derive(self, variable):
        raise NotImplementedError

class Variable(Expression):
    def __init__(self, value):
        self.value = value

    def evaluate_and_derive(self, variable):
        if self is variable:
            partial = 1.0
        else:
            partial = 0.0
        return ValueAndPartial(self.value, partial)

class Plus(Expression):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def evaluate_and_derive(self, variable):
        va = self.a.evaluate_and_derive(variable)
        vb = self.b.evaluate_and_derive(variable)
        value = va.value + vb.value
        partial = va.partial + vb.partial
        return ValueAndPartial(value, partial)

class Multiply(Expression):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def evaluate_and_derive(self, variable):
        va = self.a.evaluate_and_derive(variable)
        vb = self.b.evaluate_and_derive(variable)
        value = va.value * vb.value
        partial = vb.value * va.partial + va.value * vb.partial
        return ValueAndPartial(value, partial)

# Example: Finding the partials of z = x * (x + y) + y * y at (x, y) = (2, 3)
x = Variable(2)
y = Variable(3)
z = Plus( Multiply(x, Plus(x,y)), Multiply(y, y) )

x_partial = z.evaluate_and_derive(x).partial
y_partial = z.evaluate_and_derive(y).partial

logger.info(f"(auto) ∂z/∂x = {x_partial}, ∂z/∂y = {y_partial}")

# Example implemented using numerical differentiation of the Toolbox.
def f(_x,_y):
    return _x*(_x+_y) + _y*_y
def fy(_x):
    return f(_x,3)
def fx(_y):
    return f(2,_y)

x_partial = diff_deriv(2,fy)
y_partial = diff_deriv(3,fx)

logger.info(f" (num) ∂z/∂x = {x_partial}, ∂z/∂y = {y_partial}")
