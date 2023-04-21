import numpy as np
import math
from scipy.integrate import quad, trapz
from typing import Callable
from scipy.constants import h, c, k
from math import pi
import matplotlib.pyplot as plt
from calculating import *


def U_F(ε_values, τ_values, s_values, M_λ_values, λ1, λ2, q, D_in, f, R_l, R_i, K_dai, S_λm, T):
    # Define the integrand function
    integrand = s_values * ε_values * τ_values * M_λ_values
    
    # Compute the flux
    flux = np.trapz(integrand, dx=1e-9)
    
    # Compute the voltage at the output of the photocurrent amplifier
    U_F = (q / 4) * (D_in / f)**2 * ((S_λm * R_l * R_i * K_dai) / (R_l + R_i)) * flux
    
    return U_F

# Define the input arrays
ε_values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
τ_values = np.array([0.8, 0.85, 0.9, 0.95, 1.0])
s_values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
M_λ_values = np.array([1e-6, 2e-6, 3e-6, 4e-6, 5e-6])
λ1 = 400e-9
λ2 = 700e-9
q = 0.01
D_in = 0.1
f = 1
R_l = 1000
R_i = 100
K_dai = 10
S_λm = 1

# Define the temperature
T = 1000 # in units of K

# Calculate U_F
U_F = U_F(ε_values, τ_values, s_values, M_λ_values, λ1, λ2, q, D_in, f, R_l, R_i, K_dai, S_λm, T)
print(U_F)