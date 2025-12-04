import numpy as np
from scipy.optimize import root_scalar
import tristanVis.isolde as isolde
import warnings
import matplotlib.pyplot as plt

# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     fxn()




import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Base input file skeleton")
args = parser.parse_args()
input_file_name = args.input

# out_dir = "../vault/output_psi0.5_mul10_20250203_121432/"
# out_dir = "../slurm_scripts/"
# """
# Typical stuff read from the input file
# """
# input_file_name = out_dir + "temp_input_psi0.5_mul10.in" 

input_params    = isolde.parseInput(input_file_name)
m1              = input_params["particles"]["m1"]
m2              = input_params["particles"]["m2"]
interval        = input_params["output"]["interval"]
lst_time        = input_params["time"]["last"]
grid_x          = int(input_params["grid"]["mx0"])
grid_y          = int(input_params["grid"]["my0"])
Nsteps          = int(lst_time// interval)
CC              = input_params["algorithm"]["c"]
COMP            = input_params["plasma"]["c_omp"]
ppc0            = input_params["plasma"]["ppc0"]
SIGMA           = input_params["plasma"]["sigma"]

# weight_fac      = 2.405980
# omegap0         = CC / COMP * np.sqrt(weight_fac)
# time            = hist["time"] * omegap0
B_0             = input_params["problem"]["B_0"]#["B_0""B_field"]
B_amp           = input_params["problem"]["B_amplitude"]
duration        = input_params["problem"]["duration"]
freq            = input_params["problem"]["frequency"]
# spread          = input_params["problem"]["spread"]
# mode            = input_params["problem"]["mode"]
psi             = input_params["problem"]["psi"]
mult1           = input_params["problem"]["multiplicity_1"]
mult2           = input_params["problem"]["multiplicity_2"]
wA_wp           = input_params["problem"]["wA_wp"]
ramp_width      = input_params["problem"]["ramp_width"]
wall_x_location = input_params["problem"]["wall_x"]

mode = freq / (CC * np.cos(psi))
init_x_boundary = 5 * np.pi / mode + wall_x_location
fin_x_boundary = init_x_boundary + ramp_width * 2 * np.pi / mode



# print(wA_wp * 2 * np.pi / mode)

# def solve_for_beta(k):
#     def f(x):
#         y = 2 - x
#         return (1.0 / x**2) - (1.0 / y**2) - 2.0 * k
#     x = root_scalar(f, bracket = [1e-15, 2 - 1e-15], x0 = 1).root
#     y = 2 - x
#     beta_p = (1 - x**2)/(1 + x**2)
#     beta_m = (1 - y**2)/(1 + y**2)
#     return beta_p, beta_m

# def b_wave(x):
#     if x <= 4 * np.pi / mode:
#         return (B_amp/0.869619) * np.sin(x * mode) * np.sin(x * mode/4)**2
#     return 0

# def db_dx_wave(x):
#     if x <= 4 * np.pi / mode:
#         return -(B_amp/0.869619) * 0.125*(mode*(np.cos((mode*x)/2)- 4*np.cos(mode*x) + 3*np.cos((3*mode*x)/2)))
#     return 0

# def bg_density_scaled(x):
#     n_0_local = mult1
#     if mult1 != mult2:        
#         k = np.log((mult1-mult2) / (0.01 * mult2))/ (fin_x_boundary - init_x_boundary)
#         if x < init_x_boundary:
#             n_0_local = mult1
#         else:
#             n_0_local = mult1 + (mult2 - mult1) * (1.0 - (1.0 + k*(x - init_x_boundary)) * np.exp(-k*(x - init_x_boundary)))
#         # n_0_local = mult1 * np.exp(- k * (x - init_x_boundary))
#     # if (x <= init_x_boundary):
#     #     n_0_local = mult1
    
#     # elif (x > fin_x_boundary):
#     #     n_0_local = mult2

#     # else :
#     #     t = (x - init_x_boundary)/(fin_x_boundary - init_x_boundary)        
#     #     n_0_local = mult1 + (mult2 - mult1) * (3.0*t**2 - 2.0*t**3)
    
    
#     return n_0_local

def s_depth(COMP, show_det = False):
    B_norm          = CC**2 * SIGMA**0.5 / COMP
    unit_ch         = CC**2 / (ppc0 * COMP**2)
    weight = (mode * np.sin(psi) * B_norm/unit_ch) / (0.5 * ppc0)
    skin_depth = COMP / np.sqrt(weight)
    if show_det:
        # print("jA_ec_max = ", jA_ec_max)
        print("weight = ", weight)
        print("skin_depth = ", skin_depth)
        
    return skin_depth

# print("target s_depth = ", wA_wp * 2 * np.pi / mode)

def opt_sd(COMP):
    return s_depth(COMP) - wA_wp * 2 * np.pi / mode


with warnings.catch_warnings(action="ignore"):
    
    COMP_true = root_scalar(opt_sd, method = "brentq", bracket=[1e-3, 10 * COMP]).root

    
    # print("optimised comp = ", opt_sd(COMP_true))
    # s_depth(COMP_true, True)
    # comps = np.linspace(1e-3, 100, 1001)
    # plt.plot(comps, opt_sd(comps))
    # plt.show()

    print(COMP_true)
    # s_depth(COMP_true, True)
