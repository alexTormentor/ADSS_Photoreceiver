import numpy as np
import math
from scipy.integrate import quad
from typing import Callable
from scipy.constants import h, c, k
from math import pi
import matplotlib.pyplot as plt

def blackbody_radiation(lam, T):
    c2 = 1.4388e-2 # m K
    x_lam = 4.9651 * (lam * T / c2)
    return 142.34 * math.pow(x_lam, -5) * math.pow((math.exp(4.9651 / x_lam) - 1), -1)

def F(q, D_in, f, τ, M_λ, ε, λ1, λ2, T):
    # Define the integrand function
    integrand = lambda x: ε(x, T) * τ(x) * M_λ(x, T)
    
    # Compute the flux
    flux, error = quad(integrand, λ1, λ2)

    return (q / 4) * (D_in / f)**2 * flux

def U_F(q, D_in, f, S_λm, R_l, R_i, K_dai, integrand, λ1, λ2):
    """
    Calculates the voltage at the output of the photocurrent amplifier.
    
    Args:
    q (float): area of the slit aperture
    D_in (float): diameter of the entrance pupil
    f (float): focal length of the lens
    S_λm (float): maximum spectral sensitivity of the photodetector
    R_l (float): load resistance
    R_i (float): internal resistance of the photodetector
    K_dai (float): gain of the photocurrent amplifier
    integrand (function): function representing the integrand of the formula
    λ1 (float): lower bound of the wavelength range of the integral
    λ2 (float): upper bound of the wavelength range of the integral
    
    Returns:
    float: the voltage at the output of the photocurrent amplifier
    """
    S_λ = lambda x: 0.5
    τ = lambda x: 0.9
    M_λ = lambda x, y: 1e-6
    ε = lambda x, y: 0.5
    s = lambda λ: S_λ(λ) / S_λm
    
    # Compute the integral of the integrand function
    flux, error = quad(integrand, λ1, λ2)
    
    U_F = (q / 4) * (D_in / f)**2 * ((S_λm * R_l * R_i * K_dai) / (R_l + R_i)) * flux
    
    return U_F