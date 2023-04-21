import numpy as np
import math
from scipy.integrate import quad, trapz
from typing import Callable
from scipy.constants import h, c, k, e
from math import pi
import matplotlib.pyplot as plt
from calculating import *
import tests


def photoresistor_total_current(V_n, l, d, σ_0, Δσ, dH):
    """
    Calculates the total current through a photoresistor.
    
    Parameters:
    - V_n: supply voltage (in Volts)
    - l: length of the photoresistor (in meters)
    - d: width of the photoresistor (in meters)
    - σ_0: conductivity at zero light intensity (in Siemens per meter)
    - Δσ: change in conductivity due to light (in Siemens per meter)
    - dH: thickness of the photosensitive element (in meters)
    
    Returns:
    - The total current through the photoresistor (in Amperes)
    """
    # Calculate the two terms in the total current formula
    term1 = σ_0 * (dH / l) * V_n
    term2 = Δσ * (dH / l) * V_n
    
    # Add the two terms to get the total current
    I = term1 + term2
    
    return I

def photocurrent(e, h, c, alpha, n, tau, ue, uh, F, Vn, H, l, wavelength):
    """
    Calculates the photocurrent in a semiconductor photodetector.
    
    Parameters:
    - e: electron charge in Coulombs
    - h: thickness of the semiconductor wafer in meters
    - c: propagation velocity of electromagnetic oscillations in meters per second
    - alpha: spectral absorption index in meters^-1
    - n: quantum output of the semiconductor in photons^-1
    - tau: lifetime of carriers in seconds
    - ue: mobility of electrons in meters^2 per Volt-second
    - uh: mobility of holes in meters^2 per Volt-second
    - F: spectral density of the radiation flux in Watts per meter^-2 per meter
    - Vn: supply voltage in Volts
    - H: width of the semiconductor wafer in meters
    - l: length of the semiconductor wafer in meters
    
    Returns:
    - The photocurrent in the semiconductor photodetector (in Amperes)
    """
    k = 1.38e-23   # Boltzmann constant in Joules per Kelvin
    T = 300        # Temperature in Kelvin
    
    # Compute the spectral sensitivity
    spectral_sensitivity = alpha * n * wavelength * (e / (h * c)) * tau * (ue + uh)
    
    # Compute the flux
    flux = (1 / h) * F
    
    # Compute the photocurrent
    photocurrent = spectral_sensitivity * flux * Vn * (H / l)**2
    
    return photocurrent

def dark_current(V, T, Isat, A):
    """
    Calculates the dark current in a photodetector.
    
    Parameters:
    - V: voltage across the photodetector (in volts)
    - T: temperature (in Kelvin)
    - Isat: saturation current (in Amperes)
    - A: coefficient that depends on the material of the photodetector (default value is 1)
    
    Returns:
    - The dark current in the photodetector (in Amperes)
    """
    q = 1.6e-19   # Elementary charge in Coulombs
    k = 1.38e-23  # Boltzmann constant in Joules per Kelvin
    
    return Isat * (math.exp((-q * V) / (k * A * T)) - 1)

def current_sensitivity(e, h, c, alpha, n, tau, ue, uh, Vn, H, l, wavelength):
    """
    Calculates the monochromatic current sensitivity of a photodetector.
    
    Parameters:
    - e: electron charge in Coulombs
    - h: thickness of the semiconductor wafer in meters
    - c: propagation velocity of electromagnetic oscillations in meters per second
    - alpha: spectral absorption index in meters^-1
    - n: quantum output of the semiconductor in photons^-1
    - tau: lifetime of carriers in seconds
    - ue: mobility of electrons in meters^2 per Volt-second
    - uh: mobility of holes in meters^2 per Volt-second
    - Vn: supply voltage in Volts
    - H: width of the semiconductor wafer in meters
    - l: length of the semiconductor wafer in meters
    - wavelength: wavelength of the incident light in meters
    
    Returns:
    - The monochromatic current sensitivity of the photodetector (in Amperes per Watt)
    """
    k = 1.38e-23   # Boltzmann constant in Joules per Kelvin
    T = 300        # Temperature in Kelvin
    
    # Compute the spectral sensitivity
    spectral_sensitivity = alpha * n * wavelength * (e / (h * c)) * tau * (ue + uh)
    
    # Compute the monochromatic current sensitivity
    current_sensitivity = spectral_sensitivity * Vn * (H / l)**2
    
    return current_sensitivity

# Define the input parameters
V_n = 1.0     # Volts
l = 0.1       # meters
d = 0.01      # meters
σ_0 = 0.01    # Siemens/meter
Δσ = 0.1      # Siemens/meter
dH = 2e-3     # meters

# Calculate the total current through the photoresistor
I = photoresistor_total_current(V_n, l, d, σ_0, Δσ, dH)

# Print the result
print("Total current through the photoresistor:", I, "A")

# Define the parameters
h = 1e-6       # Thickness of the semiconductor wafer in meters
alpha = 1e5    # Spectral absorption index in meters^-1
n = 1e20       # Quantum output of the semiconductor in photons^-1
tau = 1e-7     # Lifetime of carriers in seconds
ue = 1000      # Mobility of electrons in meters^2 per Volt-second
uh = 500       # Mobility of holes in meters^2 per Volt-second
F = 1e-6       # Spectral density of the radiation flux in Watts per meter^-2 per meter
Vn = 5         # Supply voltage in Volts
H = 2e-3       # Width of the semiconductor wafer in meters
l = 10e-3      # Length of the semiconductor wafer in meters
wavelength = 5e-6

# Call the photocurrent function to compute the photocurrent
I_photocurrent = photocurrent(e, h, c, alpha, n, tau, ue, uh, F, Vn, H, l, wavelength)
print("Photocurrent:", I_photocurrent, "A")

# Calculate the dark current
I_dark = dark_current(0.1, 300, 1e-12, 4)*-1
print("Dark current: ", I_dark)

# Add the dark current and photocurrent together
I_total = I_dark + I_photocurrent
print("Total current: ", I_total)

Current_Sens = current_sensitivity(e, h, c, alpha, n, tau, ue, uh, Vn, H, l, wavelength)

print("monochromatic current sensitivity of a photodetector = ", I_photocurrent / F)

#print("monochromatic current sensitivity of a photodetector = ", I_photocurrent / Current_Sens)