latqcdtools.testing
=============

`concludeTest(lpass)`


`gaudif_results(res, res_err, res_true, res_err_true, text='', qcut=0.05)`

Compares element-by-element the results of res with res_true using Gaussian difference test, i.e. it checks
to see whether res and res_true are statistically compatible. 

`print_results(res, res_true, res_err=None, res_err_true=None, text='', prec=1e-10, abs_prec=None)`

Compares element-by-element the results of res with res_true. (Does the same with res_err and res_err_true,
if you like.) Carries out with precision prec. Use abs_prec for comparisons with zero. 

`print_results_iter(res, res_true, text)`


