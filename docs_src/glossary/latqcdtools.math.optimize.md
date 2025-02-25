latqcdtools.math.optimize
=============

`minimize(func, jac=None, hess=None, start_params=None, tol=1e-12, maxiter=10000, algorithm=None)`

Wrapper for scipy.optimize.minimize. Helps that all algorithms have common syntax.

Args:
    func: to-be-minimized function 
    jac (func, optional): Explicit Jacobian. Defaults to None.
    hess (func, optional): Explicit Hessian. Defaults to None.
    start_params (array-like, optional): Starting guess. Defaults to None.
    tol (float, optional): Solve tolerance. Defaults to 1e-12.
    maxiter (int, optional): Maximum number of solve iterations. Defaults to 10000.
    algorithm (str, optional): Solve algorithm. Defaults to Scipy default, usually BFGS.

Raises:
    e: scipy.optimize exception 
    ValueError: minimizer fails to converge

Returns:
    array-like: solution vector 

`persistentSolve(LHS, guess, tol=1e-08, maxiter=300)`

Attempt to solve LHS==0 using, in this order, SciPy's newton_krylov, fsolve, and root.

Args:
    LHS (func)
    guess: initial guess for solution 
    tol (real, optional): Solve tolerance. Defaults to 1e-8.
    maxiter (int, optional): Maximum iterations. Defaults to 300.

`solve(LHS, guess, tol=1e-08, maxiter=300, method='newton_krylov')`

Wrapper for various methods to solve LHS==0. This is to simplify the interface.

Args:
    LHS (func)
    guess: initial guess for solution 
    tol (real, optional): Solve tolerance. Defaults to 1e-8.
    maxiter (int, optional): Maximum iterations. Defaults to 300.
    method (str, optional): Defaults to 'newton_krylov'.

