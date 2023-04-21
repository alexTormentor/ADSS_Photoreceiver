import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import math

# Constants
h = 6.626e-34  # Planck's constant
c = 3e8  # Speed of light
k = 1.38e-23  # Boltzmann constant

# Functions
pi = np.pi

def planck(wav, T):
    a = 2.0*h*pi*c**2
    b = h*c/(wav*k*T)
    intensity = a/ ( (wav**5)*(math.e**b - 1.0) )
    return intensity


def tau(lambdas):
    """
    Transmission coefficient for the silicate material
    """
    return 0.95 * np.exp(-((lambdas - 2.5) / 0.1)**2)

def chi_Si(lambdas):
    """
    Spectral efficiency of the silicate material
    """
    y1_lambda = tau(lambdas) * planck(lambdas, 873)
    y2_lambda = tau(lambdas) * planck(lambdas, 1273)
    integrand_num = y1_lambda**2
    integrand_denom = y1_lambda
    numerator = simps(integrand_num, lambdas)
    denominator = simps(integrand_denom, lambdas)
    return integrand_num * numerator / denominator

def chi_FSA(lambdas):
    """
    Spectral efficiency of the FSA photoresistor
    """
    S_lambda = np.ones_like(lambdas)
    S_lambda[(lambdas < 0.4) | (lambdas > 1.1)] = 0
    return S_lambda

def chi_InSb(lambdas):
    """
    Spectral efficiency of InSb
    """
    S_lambda = np.zeros_like(lambdas)
    S_lambda[(lambdas > 1.7) & (lambdas < 5)] = 0.7
    return S_lambda

def chi_Ge(lambdas):
    """
    Spectral efficiency of Ge single crystal:Hg at 30K
    """
    S_lambda = np.zeros_like(lambdas)
    S_lambda[(lambdas > 1.8) & (lambdas < 4)] = 0.6
    return S_lambda

# Define wavelengths
lambdas = np.linspace(0.5 * 2.28, 3 * 3.32, 1000)

# Calculate spectral efficiencies
chi_Si_lambda = chi_Si(lambdas)
chi_FSA_lambda = chi_FSA(lambdas)
chi_InSb_lambda = chi_InSb(lambdas)
chi_Ge_lambda = chi_Ge(lambdas)

# Plot results
fig, ax = plt.subplots()
ax.plot(lambdas, chi_Si_lambda, label='Silicate material')
ax.plot(lambdas, chi_FSA_lambda, label='FSA photoresistor')
ax.plot(lambdas, chi_InSb_lambda, label='InSb')
ax.plot(lambdas, chi_Ge_lambda, label='Ge single crystal:Hg at 30K')
ax.set_xlabel('Wavelength (μm)')
ax.set_ylabel('Spectral efficiency')
ax.legend()
plt.show()


def tau_FSA(lambdas):
    """
    Transmission coefficient for the FSA photoresistor
    """
    return 0.6 * np.exp(-((lambdas - 0.9) / 0.1)**2)

def tau_InSb(lambdas):
    """
    Transmission coefficient for InSb
    """
    return 0.9 * np.exp(-((lambdas - 3.6) / 0.1)**2)

def tau_Ge(lambdas):
    """
    Transmission coefficient for Ge single crystal:Hg at 30K
    """
    return 0.8 * np.exp(-((lambdas - 2.5) / 0.1)**2)

def tau_Psb(lambdas):
    """
    Transmission coefficient for Psb material
    """
    return 0.8 * np.exp(-((lambdas - 1.5) / 0.05)**2)

# Define wavelengths
lambdas = np.linspace(0.5 * 2.28, 3 * 3.32, 1000)

# Calculate transmission coefficient
tau_lambda = tau(lambdas)
tau_lambda_FSA = tau_FSA(lambdas)
tau_lambda_InSb = tau_InSb(lambdas)
tau_lambda_Ge = tau_Ge(lambdas)
tau_lambda_Psb = tau_Psb(lambdas)

# Plot the transmission coefficient
fig, ax = plt.subplots()
#ax.plot(lambdas, tau_lambda, label='Silicate material')
#ax.plot(lambdas, tau_lambda_FSA, label='FSA photoresistor')
ax.plot(lambdas, tau_lambda_InSb, label='InSb')
ax.plot(lambdas, tau_lambda_Ge, label='Ge single crystal:Hg at 30K')
ax.plot(lambdas, tau_lambda_Psb, label='Psb')
ax.set_xlabel('Wavelength (microns)')
ax.set_ylabel('Transmission coefficient')
ax.legend()
plt.show()


def s_Si(lambdas):
    """
    Spectral efficiency of the silicate material
    """
    y1_lambda = tau(lambdas) * planck(lambdas, 873)
    y2_lambda = tau(lambdas) * planck(lambdas, 1273)
    integrand_num = y1_lambda**2
    integrand_denom = y2_lambda
    numerator = simps(integrand_num, lambdas)
    denominator = simps(integrand_denom, lambdas)
    return integrand_num * numerator / denominator

def s_FSA(lambdas):
    """
    Spectral efficiency of the FSA photoresistor
    """
    S_lambda = np.zeros_like(lambdas)
    S_lambda[(lambdas > 0.4) & (lambdas < 1.1)] = 0.5
    return S_lambda * tau_FSA(lambdas)

def s_InSb(lambdas):
    """
    Spectral efficiency of InSb
    """
    S_lambda = np.zeros_like(lambdas)
    S_lambda[(lambdas > 1.7) & (lambdas < 5)] = 0.7
    return S_lambda * tau_InSb(lambdas)

def s_Ge(lambdas):
    """
    Spectral efficiency of Ge single crystal:Hg at 30K
    """
    S_lambda = np.zeros_like(lambdas)
    S_lambda[(lambdas > 1.8) & (lambdas < 4)] = 0.6
    return S_lambda * tau_Ge(lambdas)

def s_Psb(lambdas):
    """
    Spectral efficiency of Psb at 295K
    """
    S_lambda = np.zeros_like(lambdas)
    S_lambda[(lambdas > 0.4) & (lambdas < 1.1)] = 0.4
    return S_lambda * tau_Psb(lambdas)

# Calculate the spectral efficiency for each material
s_Ge_lambda = s_Ge(lambdas)
s_InSb_lambda = s_InSb(lambdas)
s_FSA_lambda = s_FSA(lambdas)
s_Psb_lambda = s_Psb(lambdas)
S = s_Si(lambdas)

# Plot the spectral efficiency functions
plt.plot(lambdas, s_Ge_lambda, label='Ge:Hg at 30K')
plt.plot(lambdas, s_InSb_lambda, label='InSb at 295K')
#plt.plot(lambdas, s_FSA_lambda, label='FSA photoresistor')
plt.plot(lambdas, s_Psb_lambda, label='Psb at 295K')
#plt.plot(lambdas, S, label='Silicate')
plt.legend()
plt.xlabel('Wavelength (μm)')
plt.ylabel('Spectral efficiency')
plt.show()


def y1(lambdas, T):
    """
    Dimensionless value of the spectral density function of the energy luminosity at the output of the optical element
    """
    return (2 * h * c**2 / lambdas**5) * 1 / (np.exp(h * c / (lambdas * k * T)) - 1)

def y2(lambdas, T):
    """
    Dimensionless value of the spectral density function of the photon luminosity at the output of the optical element
    """
    return (2 * h * c**2 / lambdas**5) * (h * c / (lambdas * k * T)) / (np.exp(h * c / (lambdas * k * T)) - 1)

# Define a range of temperatures
temperatures = np.arange(1000, 5001, 1000)

# Define a range of wavelengths
lambdas = np.linspace(0.5 * 2.28, 3 * 3.32, 1000)

# Calculate y1 values for each temperature and wavelength
y1_values = []
for T in temperatures:
    y1_values.append(y1(lambdas, T))

# Plot y1 values for each temperature as a separate line
fig, ax = plt.subplots()
for i, T in enumerate(temperatures):
    ax.plot(lambdas, y1_values[i], label=f'T = {T} K')
ax.set_xlabel('Wavelength (μm)')
ax.set_ylabel('Dimensionless spectral density function')
ax.set_title('y1(λ,T) vs. λ for different temperatures')
ax.legend()
plt.show()

def trapezoidal(x, y, z, c, s):
    """
    Trapezoidal membership function
    """
    result = np.zeros_like(x)
    mask1 = np.logical_and(y <= x, x < c - 0.5)
    mask2 = np.logical_and(c + 0.5 < x, x <= z)
    mask3 = np.logical_and(c - 0.5 <= x, x <= c + 0.5)
    
    if np.any(y >= c + 0.5):
        result[:] = 0
    else:
        result[mask1] = 0
        result[mask2] = 0
        if np.any(y == c - 0.5):
            result[mask3] = s
        else:
            result[mask3] = s * (x[mask3] - y) / (c - y)
        r = np.any(((c - 0.5) <= x <= (c + 0.5)), None)
        result[r] = 1

    return result


# Define approximate spectral characteristics of photodetectors and optical materials
lambdas = np.linspace(0.3, 12, 1000)
tau_Psb_ = 0.9  # Example value
s_InSb_ = 1.0  # Example value
s_Ge_ = 0.6  # Example value

# Define the AFT temperatures
T_1 = 873  # 600 °C
T_2 = 1273  # 1000 °C

# Calculate the spectral density function of the energy luminosity at the output of the optical element
y_873 = trapezoidal(s_InSb_ * tau_Psb_ * s_Ge_ * planck(lambdas, T_1), s_InSb_ * tau_Psb_ * s_InSb_ * planck(lambdas, T_1), 1, 0.5, 1)
y_1273 = trapezoidal(s_InSb_ * tau_Psb_ * s_Ge_ * planck(lambdas, T_2), s_InSb_ * tau_Psb_ * s_InSb_ * planck(lambdas, T_2), 1, 0.5, 1)

# Plot the results
plt.plot(lambdas, y_873, label='873K')
plt.plot(lambdas, y_1273, label='1273K')
plt.legend()
plt.xlabel('Wavelength (µm)')
plt.ylabel('y(λ,T)')
plt.show()