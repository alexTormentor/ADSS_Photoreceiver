import numpy as np
import math
from scipy.integrate import quad, trapz
from typing import Callable
from scipy.constants import h, c, k
from math import pi
import matplotlib.pyplot as plt
from calculating import *


# Define the F function
def F(q, D_in, f, τ, M_λ, ε, λ1, λ2, T):
    # Define the integrand function
    integrand = lambda x: ε(x, T) * τ(x) * M_λ(x, T)
    
    # Compute the flux
    flux, error = quad(integrand, λ1, λ2)

    return (q / 4) * (D_in / f)**2 * flux

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

# Select the index of the value you want to use for each coefficient
index = 2

# Define the functions for each coefficient using the selected index
ε = lambda x, T: ε_values[index]
τ = lambda x: τ_values[index]
M_λ = lambda x, T: M_λ_values[index]

# Compute the output voltage U_F
U_F = F(q, D_in, f, τ, M_λ, ε, λ1, λ2, T=300)
#print(U_F)

# темновой ток
I_photo = U_F * S_λm * 0.05

I_dark = 0.002 * 1 - I_photo
print(I_dark)