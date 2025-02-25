latqcdtools.base.speedify
=============

`compile(func)`


`numbaList(inList)`

Turn a list into List that numba can parse. 

`numbaOFF()`

Turn off numba compilation for small functions. 

`numbaON()`

Use numba wherever possible. By default it is turned off, since compilation takes some time,
and hence you will only see a performance boost for particularly long-running functions. Must be
called at the beginning of your code. 

`parallel_function_eval(function, input_array, args=(), nproc=6, parallelizer='pathos.pools')`

Parallelize a function over an input_array. Effectively this can replace a loop over an array and should
lead to a performance boost.

Args:
    function (func): to-be-parallelized function 
    input_array (array-like): array over which it should run 
    nproc (int): number of processes 

Returns:
    array-like: func(input_array)

`parallel_reduce(function, input_array, args=(), nproc=6, parallelizer='pathos.pools') -> float`

Parallelize a function over an input_array, then sum over the input_array elements. 

Args:
    function (func): to-be-parallelized function 
    input_array (array-like): array over which it should run 
    nproc (int): number of processes 

Returns:
    float

`ComputationClass(function, input_array, args, nproc, parallelizer)`


