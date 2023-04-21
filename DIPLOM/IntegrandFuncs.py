import numpy as np
import math
from scipy.integrate import quad, trapz
from typing import Callable
from scipy.constants import h, c, k, e
from math import pi
import matplotlib.pyplot as plt




def spectral_efficiency(F_e_lambda, s_lambda, lambda_vals):
    """
    Calculates the spectral efficiency of a converter based on the spectral density of the input radiation flux and 
    the spectral sensitivity of the converter.
    
    Parameters:
    - F_e_lambda: spectral distribution of the spectral density of the input radiation flux
    - s_lambda: spectral sensitivity of the converter
    - lambda_vals: array of wavelength values in meters
    
    Returns:
    - The spectral efficiency of the converter
    """
    phi_e_lambda = F_e_lambda / np.max(F_e_lambda)
    integrand_1 = s_lambda * phi_e_lambda
    integrand_2 = phi_e_lambda
    
    numerator = np.trapz(integrand_1, lambda_vals)
    denominator = np.trapz(integrand_2, lambda_vals)
    
    spectral_efficiency = numerator / denominator
    
    return spectral_efficiency


# define the spectral distribution of the input radiation flux F_e(λ)
lambdas = np.linspace(300, 1100, 801)  # wavelengths in nm
F_e = 1e-3 * np.exp(-(lambdas - 700)**2 / (2*50**2))  # spectral density in W/m^2/nm

# calculate the spectral efficiency
s = np.ones_like(lambdas)  # assume perfect spectral sensitivity
chi = spectral_efficiency(s, F_e, lambdas)

print(f"Spectral efficiency: {chi:.2e}")


def spectral_efficiency2(lens_transmission, detector_sensitivity):
    """
    Calculates the spectral efficiency of a photodetector and lens system based on the transmission coefficient
    of the lens and sensitivity of the detector.

    Args:
    - lens_transmission: array-like, transmission coefficient of the lens for each wavelength
    - detector_sensitivity: array-like, sensitivity of the detector for each wavelength

    Returns:
    - spectral efficiency: float, the spectral efficiency of the system
    """

    y_1 = np.multiply(lens_transmission, detector_sensitivity)
    integral_s = np.trapz(np.multiply(detector_sensitivity, y_1))
    integral = np.trapz(y_1)

    return integral_s / integral

# define arrays for lens transmission and detector sensitivity
lens_transmission = np.array([0.8, 0.7, 0.6])
detector_sensitivity = np.array([0.5, 0.6, 0.7])

# calculate spectral efficiency
print(spectral_efficiency2(lens_transmission, detector_sensitivity))



def output_spectrum(wavelength):
    # Example output spectrum function
    return 1e-9 * wavelength ** 2

def input_spectrum(wavelength):
    # Example input spectrum function
    return 1e-6 * wavelength ** 3
def spectral_efficiency(spectrum, input_spectrum, wavelength_range):
    """
    Calculate the spectral efficiency of a converter system.

    Args:
        spectrum (callable): Function that returns the spectral density of the output radiation flux
                             at a given wavelength in units of W/m^2/nm.
        input_spectrum (callable): Function that returns the spectral density of the input radiation flux
                                   at a given wavelength in units of W/m^2/nm.
        wavelength_range (tuple): Range of wavelengths to integrate over in units of nm.

    Returns:
        The spectral efficiency of the converter system.
    """
    integrand1 = lambda x: spectrum(x) * input_spectrum(x) / input_spectrum(wavelength_range[1])
    integrand2 = lambda x: input_spectrum(x) / input_spectrum(wavelength_range[1])
    numerator, _ = quad(integrand1, wavelength_range[0], wavelength_range[1])
    denominator, _ = quad(integrand2, wavelength_range[0], wavelength_range[1])
    return numerator / denominator


wavelength_range = (600, 800)  # Example wavelength range

spectral_eff = spectral_efficiency(output_spectrum, input_spectrum, wavelength_range)
print(spectral_eff)