#!/usr/bin/env python3
import argparse
import sys
from latqcdtools.interpolate_bspline import *
from latqcdtools.plotting import *
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("beta", type=float) 

parser.add_argument("--interpolate", default = False, action="store_true") 
parser.add_argument("--wissel", default = False, action="store_true") 
parser.add_argument("--show_plot", default = False, action="store_true") 

use_luescher = True

args = parser.parse_args()

if not args.beta:
    print("Usage " + sys.argv[0] + "--beta beta")
    sys.exit(-1)

beta=args.beta

g = 6.0/beta
k1 = 0.008439857

k = 0.125 + k1*g + 0.0085*g**2 - 0.0272*g**3 + 0.0420*g**4 - 0.0204*g**5


if args.interpolate:
    print(beta, k)

betas = [6.0,
        6.2, 
        6.4, 
        6.8, 
        7.4, 
        8.0, 
        9.6]

kappa_cs = [
        0.135196,
        0.135795,
        0.135720,
        0.135097,
        0.134071,
        0.133173,
        0.131448
        ]


if args.wissel:
    betas = [6.00,
             6.136, 
             6.205, 
             6.338, 
             6.499, 
             6.640, 
             6.721,
             6.872,
             7.192,
             7.457]

    kappa_cs = [
            0.13520,
            0.13571,
            0.13580,
            0.13580,
            0.13558,
            0.13536,
            0.13522,
            0.13495,
            0.13437,
            0.13396
            ]

int_beta, int_kappa_c = spline(betas, kappa_cs, xspline=[beta])

if not args.interpolate:
    print(int_beta[0], int_kappa_c[0])

#int_beta, int_kappa_c = spline(betas, kappa_cs)

if args.show_plot:
    int_beta_cont, int_kappa_c_cont = spline(betas, kappa_cs)
    plot_dots(int_beta, int_kappa_c)
    plot_lines(int_beta_cont, int_kappa_c_cont, marker=None)
    plot_dots(betas, kappa_cs)
    plt.show()
