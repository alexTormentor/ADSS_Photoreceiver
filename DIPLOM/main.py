import numpy as np
import math
from scipy.integrate import quad
from typing import Callable
from scipy.constants import h, c, k
from math import pi
import matplotlib.pyplot as plt
from tests import *
from calculating import *

def t1():
    test3_1()
    test4_1()
    test5_1()
    test6_1()
    test7_1()
    test8_1()
    test9_1()
    test10_1()
    test11_1()
    test12_1()
    test13_1()
    test14()
    test15()
    test16()
    test17()
    test18()
    test19()
    test20()
    test21()
    test22()
    test23()
    test24()
    test25()
    test26()
    test27()
    test28()
    test30()
    test31()
    test32()
    test33()
    test34()
    test35()
    test36()
    test37()
    test38()
    test39()
    test40()
    test40_2()
    test41()
    test42()
    test43()



if __name__ == '__main__':
    #t1()

    # Test the function with λ = 500 nm and T = 1000 K
    lam = 500 # in units of nm
    T = 1000 # in units of K
    result = blackbody_radiation(lam*1e-9, T) # Convert λ from nm to m
    print(f"The spectral density function for λ = {lam} nm and T = {T} K is {result:.2e}.")

    # Test the function with λ = 1000 nm and T = 600 K
    lam = 1000 # in units of nm
    T = 600 # in units of K
    result = blackbody_radiation(lam*1e-9, T) # Convert λ from nm to m
    print(f"The spectral density function for λ = {lam} nm and T = {T} K is {result:.2e}.")
    print('\n\n\n')
    # Define values for parameters
    q = 0.01
    D_in = 0.1
    f = 1
    τ = lambda x: 0.9 # Example transmission coefficient function that returns a constant value
    M_λ = lambda x, y: 1e-6 # Example spectral density of the energy luminosity function that returns a constant value
    ε = lambda x, y: 0.5 # Example spectral emittance function that returns a constant value
    λ1 = 400e-9
    λ2 = 700e-9
    T = 1000 # in units of K

    # Compute the flux
    flux = F(q, D_in, f, τ, M_λ, ε, λ1, λ2, T)

    print(f"The flux is {flux:.2e} W/m^2.")

    print('\n\n\n')
    # Define values for parameters
    R_i = 1e4
    R_l = 1e3
    K_dai = 1e3
    S_λm = 0.8
    S_λ = lambda x: 0.5

    # Define the integrand function
    integrand = lambda x: S_λ(x) / S_λm * ε(x, T) * τ(x) * M_λ(x, T)

    # Define the temperature
    T = 1000 # in units of K

    # Compute the output voltage
    U_F = U_F(q, D_in, f, R_i, R_l, K_dai, S_λm, integrand, λ1, λ2)
    print(f"The output voltage of the photocurrent amplifier is {U_F:.2e} volts.")


