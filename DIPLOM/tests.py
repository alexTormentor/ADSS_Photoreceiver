import numpy as np
import math
from scipy.integrate import quad
from typing import Callable
from scipy.constants import h, c, k
from math import pi
import matplotlib.pyplot as plt

# Коэффициент нелинейности ВАХ
def calculate_nonlinearity_coefficient(voltage, current):
    # Convert the voltage and current data to numpy arrays
    voltage = np.array(voltage)
    current = np.array(current)
    
    # Calculate the average slope of the curve
    avg_slope = np.mean(np.diff(current) / np.diff(voltage))
    
    # Calculate the maximum deviation from the average slope
    max_deviation = np.max(np.abs((current[1:] - current[:-1]) / (voltage[1:] - voltage[:-1]) - avg_slope))
    
    # Calculate the Coefficient of Nonlinearity
    nonlinearity_coefficient = max_deviation / avg_slope
    
    return nonlinearity_coefficient
def TestCNC():
    # Define the voltage and current data points
    voltage = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    current = [0, 0.12, 0.4, 0.8, 1.5, 2.4, 3.6]

    # Calculate the Coefficient of Nonlinearity
    nonlinearity_coefficient = calculate_nonlinearity_coefficient(voltage, current)

    # Print the result
    print("Коэффициент нелинейности: ", nonlinearity_coefficient)

# Коэффициент нелинейности света
def calculate_light_nonlinearity_coefficient(intensity, current):
    # Convert the intensity and current data to numpy arrays
    intensity = np.array(intensity)
    current = np.array(current)
    
    # Calculate the average slope of the curve
    avg_slope = np.mean(np.diff(current) / np.diff(intensity))
    
    # Calculate the maximum deviation from the average slope
    max_deviation = np.max(np.abs((current[1:] - current[:-1]) / (intensity[1:] - intensity[:-1]) - avg_slope))
    
    # Calculate the Coefficient of Nonlinearity
    nonlinearity_coefficient = max_deviation / avg_slope
    
    return nonlinearity_coefficient
def testCLNC():
        # Define the light intensity and current data points
    intensity = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    current = [0, 1, 5, 13, 25, 42, 64, 90, 120, 155, 192]

    # Calculate the Coefficient of Nonlinearity for the light characteristic
    nonlinearity_coefficient = calculate_light_nonlinearity_coefficient(intensity, current)

    # Print the result
    print("Коэффициент нелинейности световой характеристики: ", nonlinearity_coefficient)

# световой поток
def calculate_radiation_flux(energy, time):
    """
    Calculates the radiation flux based on the energy of radiation and the time over which the energy is transferred.

    Parameters:
    energy (float): The energy of radiation in Joules.
    time (float): The time over which the energy is transferred in seconds.

    Returns:
    The radiation flux in Watts per square meter.
    """
    radiation_flux = energy / time
    return radiation_flux
def testCRF():
        # Define the energy of radiation in Joules and the time over which it is transferred in seconds
    energy = 1000
    time = 10

    # Calculate the radiation flux
    radiation_flux = calculate_radiation_flux(energy, time)

    # Print the result
    print("Световой поток: ", radiation_flux, "Вт/м^2")

# преобразование излучения в ток
def radiation_current_equation(V, voltage, current, intensity, energy, time):
    C = 0.25
    F = calculate_radiation_flux(energy, time)

    # Interpolate the current values at the same intensity values
    current_interp = np.interp(intensity, voltage, current)

    nonlinearity_coefficient = calculate_nonlinearity_coefficient(voltage, current)
    light_nonlinearity_coefficient = calculate_light_nonlinearity_coefficient(intensity, current_interp)
    return C * V**nonlinearity_coefficient * F**light_nonlinearity_coefficient
def RCE():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    current = [0, 0.12, 0.4, 0.8, 1.5, 2.4, 3.6]
    intensity = [0, 1, 2]
    energy = 1000
    time = 10
    Vol = 1.6

    print('ток, преобразованный из излучения: ', radiation_current_equation(Vol, voltage, current, intensity, energy, time))

# средняя чувствительность
def AverageValSensPhotor(V, voltage, current, intensity, energy, time):
    C = 0.25
    F = calculate_radiation_flux(energy, time)

    # Interpolate the current values at the same intensity values
    current_interp = np.interp(intensity, voltage, current)

    nonlinearity_coefficient = calculate_nonlinearity_coefficient(voltage, current)
    light_nonlinearity_coefficient = calculate_light_nonlinearity_coefficient(intensity, current_interp)
    return C * V**nonlinearity_coefficient * F**(light_nonlinearity_coefficient - 1)
def AVP():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    current = [0, 0.12, 0.4, 0.8, 1.5, 2.4, 3.6]
    intensity = [0, 1, 2, 3]
    energy = 100
    time = 1.5
    Vol = 0.014

    print('средняя чувствительность: ', radiation_current_equation(Vol, voltage, current, intensity, energy, time))

# дифф. чувствительность
def DeffValSensPhotor(V, voltage, current, intensity, energy, time):
    C = 0.25
    F = calculate_radiation_flux(energy, time)

    # Interpolate the current values at the same intensity values
    current_interp = np.interp(intensity, voltage, current)

    nonlinearity_coefficient = calculate_nonlinearity_coefficient(voltage, current)
    light_nonlinearity_coefficient = calculate_light_nonlinearity_coefficient(intensity, current_interp)
    return light_nonlinearity_coefficient * C * V**nonlinearity_coefficient * F**(light_nonlinearity_coefficient - 1)
def DVP():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    current = [0, 0.12, 0.4, 0.8, 1.5, 2.4, 3.6]
    intensity = [0, 1, 2, 3]
    energy = 100
    time = 3
    Vol = 0.014

    print('дифференцированная чувствительность: ', DeffValSensPhotor(Vol, voltage, current, intensity, energy, time))

# инф. сост. тока
def calculate_information_component(V, voltage, current, intensity, energy, time):
    """
    Calculates the information component of the current.

    Parameters:
    V (float): The applied voltage in volts.
    voltage (list): The list of voltage values in volts.
    current (list): The list of current values in amperes.
    intensity (list): The list of intensity values in W/m^2.
    energy (float): The energy of radiation in Joules.
    time (float): The time over which the energy is transferred in seconds.

    Returns:
    The information component of the current.
    """
    F = calculate_radiation_flux(energy, time)
    current_equation = radiation_current_equation(V, voltage, current, intensity, energy, time)
    information_component = current_equation * F
    return information_component
def test1():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    current = [0, 0.12, 0.4, 0.8, 1.5, 2.4, 3.6]
    intensity = [0, 1, 2, 3]
    energy = 500
    time = 1.5
    Vol = 0.014

    information_component = calculate_information_component(Vol, voltage, current, intensity, energy, time)
    print('Информационная составляющая тока: ', information_component)

# темновой ток
def calculate_dark_current(V, T, Is):
    """
    Calculates the dark current in a photodiode at a given voltage and temperature.

    Parameters:
    V (float): The voltage applied to the photodiode in volts (V).
    T (float): The temperature in kelvin (K).
    Is (float): The saturation current in amperes (A).

    Returns:
    The dark current in amperes (A).
    """
    # Define the constants q and k
    q = 1.60217662e-19  # Elementary charge in coulombs
    k = 1.380649e-23  # Boltzmann constant in joules per kelvin
    
    # Calculate the dark current using the formula
    Id = Is * (np.exp(q*V/(k*T)) - 1)
    
    return Id
def test2():
    # Define the input parameters
    V = 0.5  # Voltage in volts
    T = 300  # Temperature in kelvin
    Is = 1e-9  # Saturation current in amperes

    # Calculate the dark current
    Id = calculate_dark_current(V, T, Is)

    # Print the result
    print(f"Темновой ток = {Id:.3e} А")

# преобразование светового потока в эл. ток
def photoresistor_current(C, U, a, F, b):
    """
    Calculates the electric current flowing through a photoresistor based on the given parameters.
    
    Parameters:
    - C: constant specific to the photoresistor
    - U: voltage across the photoresistor (in Volts)
    - a: coefficient of nonlinearity of the volt-ampere characteristic of the photoresistor
    - F: radiation flux (in Watts per square meter)
    - b: coefficient of nonlinearity of the light characteristic of the photoresistor
    
    Returns:
    - The electric current flowing through the photoresistor (in Amperes)
    """
    return C * (U ** a) * (F ** b)
def test3_1():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5]
    current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
    intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    C = 0.1
    U = 5.0
    a = calculate_nonlinearity_coefficient(voltage, current)
    F = calculate_radiation_flux(energy=100.0, time=1.5)
    b = calculate_light_nonlinearity_coefficient(intensity, current)

    current = photoresistor_current(C, U, a, F, b)
    print("Ток, преобразованный из светового потока:", current, "A")
def test3():
    C = 0.1
    U = 5
    a = 0.5
    F = 100
    b = 0.3

    current = photoresistor_current(C, U, a, F, b)
    print("Current:", current, "A")

# среднее значение чувствительности
def photoresistor_sensitivity(C, U, a, F, b):
    """
    Calculates the average sensitivity of a photoresistor over the range of light flux from F1 to F2.
    
    Parameters:
    - C: constant specific to the photoresistor
    - U: voltage across the photoresistor (in Volts)
    - a: coefficient of nonlinearity of the volt-ampere characteristic of the photoresistor
    - F: radiation flux (in Watts per square meter)
    - b: coefficient of nonlinearity of the light characteristic of the photoresistor
    
    Returns:
    - The average sensitivity of the photoresistor (in Amperes per Watt per square meter)
    """
    return C * (U ** a) * (F ** (b - 1))
def test4_1():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5]
    current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
    intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    C = 0.1
    U = 5
    a = calculate_nonlinearity_coefficient(voltage, current)
    F = calculate_radiation_flux(energy=100.0, time=1.5)
    b = calculate_light_nonlinearity_coefficient(intensity, current)

    sensitivity = photoresistor_sensitivity(C, U, a, F, b)
    print("Средняя чувствительность:", sensitivity, "А/Вт/м^2")
def test4():
    C = 0.2
    U = 5
    a = 0.5
    F = 100
    b = 0.3

    sensitivity = photoresistor_sensitivity(C, U, a, F, b)
    print("Sensitivity:", sensitivity, "A/W/m^2")

# дифференциальное значение чувствительности
def photoresistor_differential_sensitivity(C, U, a, F, b):
    """
    Calculates the differential sensitivity of a photoresistor.
    
    Parameters:
    - C: constant specific to the photoresistor
    - U: voltage across the photoresistor (in Volts)
    - a: coefficient of nonlinearity of the volt-ampere characteristic of the photoresistor
    - F: radiation flux (in Watts per square meter)
    - b: coefficient of nonlinearity of the light characteristic of the photoresistor
    
    Returns:
    - The differential sensitivity of the photoresistor (in Amperes per Watt per square meter)
    """
    return b * C * (U ** a) * (F ** (b - 1))
def test5_1():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5]
    current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
    intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    C = 0.1
    U = 5
    a = calculate_nonlinearity_coefficient(voltage, current)
    F = calculate_radiation_flux(energy=100.0, time=1.5)
    b = calculate_light_nonlinearity_coefficient(intensity, current)

    differential_sensitivity = photoresistor_differential_sensitivity(C, U, a, F, b)
    print("Дифференцированная чувствительность:", differential_sensitivity, "А/Вт/м^2")
def test5():
    C = 0.2
    U = 5
    a = 0.5
    F = 100
    b = 0.3

    differential_sensitivity = photoresistor_differential_sensitivity(C, U, a, F, b)
    print("Differential sensitivity:", differential_sensitivity, "A/W/m^2")

# темновой ток
def photoresistor_dark_current(V, T, I_sat, A):
    """
    Calculates the dark current in a photoresistor.
    
    Parameters:
    - V: voltage across the photoresistor (in Volts)
    - T: temperature (in Kelvin)
    - I_sat: saturation current (in Amperes)
    - A: coefficient that depends on the material of the photodetector (default value is 1)
    
    Returns:
    - The dark current in the photoresistor (in Amperes)
    """
    q = 1.6e-19   # Elementary charge in Coulombs
    k = 1.38e-23  # Boltzmann constant in Joules per Kelvin
    
    return I_sat * (math.exp((-q * V) / (k * A * T)) - 1)
def test6_1():
    V = 5    # Voltage across the photoresistor in Volts
    T = 300  # Temperature in Kelvin
    I_sat = 1e-3  # Saturation current in Amperes
    A = 4   # Coefficient for the material of the photodetector

    dark_current = photoresistor_dark_current(V, T, I_sat, A)
    print("Темновой ток:", dark_current*-1, "A")
def test6():
    V = 5    # Voltage across the photoresistor in Volts
    T = 300  # Temperature in Kelvin
    I_sat = 1e-9  # Saturation current in Amperes
    A = 4   # Coefficient for the material of the photodetector

    dark_current = photoresistor_dark_current(V, T, I_sat, A)
    print("Dark current:", dark_current, "A")

# преобразование фотодиода
def photodiode_transformations(C, U, a, F, b, T, I_sat, A):
    photSens = photoresistor_sensitivity(C, U, a, F, b)
    darkCur = photoresistor_dark_current(U, T, I_sat, A)
    return (F * photSens) - darkCur
def test7_1():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5]
    current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
    intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    C = 0.1
    U = 5
    a = calculate_nonlinearity_coefficient(voltage, current)
    F = calculate_radiation_flux(energy=100.0, time=1.5)
    b = calculate_light_nonlinearity_coefficient(intensity, current)
    T = 300  # Temperature in Kelvin
    I_sat = 1e-3  # Saturation current in Amperes
    A = 4   # Coefficient for the material of the photodetector

    trans = photodiode_transformations(C, U, a, F, b, T, I_sat, A)
    print('преобразование фотодиода: ', trans)
def test7():
    C = 0.1
    U = 5
    a = 0.5
    F = 100
    b = 0.3
    T = 300  # Temperature in Kelvin
    I_sat = 1e-9  # Saturation current in Amperes
    A = 4   # Coefficient for the material of the photodetector

    trans = photodiode_transformations(C, U, a, F, b, T, I_sat, A)
    print(trans)

# зависимость напряжения холостого хода на фотодиоде
def photoresistor_voltage(F, A, T, I_sat, C, U, a, b):
    k = 1.38064852e-23 # Boltzmann constant
    q = 1.60217662e-19 # Elementary charge
    photSens = photoresistor_sensitivity(C, U, a, F, b) * F
    return ((k * A * T) / q) * math.log(((photSens) / I_sat) + 1)
def test8_1():
    voltage = [0, 0.5, 1, 1.5, 2, 2.5]
    current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
    intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    C = 0.1
    U = 5
    a = calculate_nonlinearity_coefficient(voltage, current)
    F = calculate_radiation_flux(energy=100.0, time=1.5)
    b = calculate_light_nonlinearity_coefficient(intensity, current)
    T = 300  # Temperature in Kelvin
    I_sat = 1e-3  # Saturation current in Amperes
    A = 4   # Coefficient for the material of the photodetector

    photVolt = photoresistor_voltage(F, A, T, I_sat, C, U, a, b)
    print('зависимость напряжения холостого хода на фотодиоде: ', photVolt)
def test8():
    C = 0.1
    U = 5
    a = 0.5
    F = 100
    b = 0.3
    T = 300  # Temperature in Kelvin
    I_sat = 1e-9  # Saturation current in Amperes
    A = 4   # Coefficient for the material of the photodetector

    photVolt = photoresistor_voltage(F, A, T, I_sat, C, U, a, b)
    print(photVolt)

# фотодиодный режим
def calculate_max_product(S, Sl, Un, Is):
    U_max = (Sl / S) * (Un / Is)
    S_U_max = S * U_max
    return S_U_max
def test9_1():
    S = 0.8  # Amperes per Watt
    Sl = 0.9  # Amperes per Watt
    Un = 0.5  # Volts
    Isat = 1e-3  # Amperes
    U_max = calculate_max_product(S, Sl, Un, Isat)
    print('Максимальное напряжение в фотодиодном режиме: ', U_max)
def test9():
    S = 0.8  # Amperes per Watt
    Sl = 0.9  # Amperes per Watt
    Un = 0.5  # Volts
    Isat = 0.00001  # Amperes
    U_max = calculate_max_product(S, Sl, Un, Isat)
    print(U_max)

# фотогальванический режим
def calc_S_Umax(S, R_0, S_l, I_f, I_s):
    """
    Calculates the product of sensitivity and maximum output voltage in photovoltaic mode of a photodiode.

    Args:
    S (float): average sensitivity of the photodiode (in Amps per Watt)
    R_0 (float): load resistance of the circuit (in Ohms)
    S_l (float): spectral sensitivity of the photodiode (in Amps per Watt)
    I_f (float): photocurrent generated by the photodiode under illumination (in Amperes)
    I_s (float): dark current generated by the photodiode in the absence of illumination (in Amperes)

    Returns:
    float: the product of sensitivity and maximum output voltage (in Volts)

    """
    U_max = (R_0 * S_l) * (1 / ((I_f / I_s) + 1))
    return S * U_max
def test10_1():
    # Example input values
    S = 0.8  # Amps per Watt
    R_0 = 1000  # Ohms
    S_l = 0.9  # Amps per Watt
    I_f = 0.001  # Amperes
    I_s = photoresistor_dark_current(V=5, T=300, I_sat=1e-3, A=4)*-1  # Amperes

    # Call the function
    result = calc_S_Umax(S, R_0, S_l, I_f, I_s)

    # Print the result
    print("В фотогальваническом режиме = {:.4f} Вольт".format(result))
def test10():
    # Example input values
    S = 0.8  # Amps per Watt
    R_0 = 1000  # Ohms
    S_l = 0.9  # Amps per Watt
    I_f = 0.001  # Amperes
    I_s = 0.00001  # Amperes

    # Call the function
    result = calc_S_Umax(S, R_0, S_l, I_f, I_s)

    # Print the result
    print("S * U_max = {:.4f} Volts".format(result))

# зависимость изм. коэф. усиления от изм. напр. питания(умножители)
def delta_gain(n, delta_U_n, U_n, M):
    """
    Calculates the change in gain for a photomultiplier with n cascades given a change in voltage of the nth stage.
    
    Args:
    - n (int): the number of cascades
    - delta_U_n (float): the change in voltage of the nth stage
    - U_n (float): the original voltage of the nth stage
    - M (float): the original gain of the photomultiplier
    
    Returns:
    - delta_M (float): the change in gain
    """
    delta_M = (0.7 / 10) * n * (delta_U_n / U_n) * M
    return delta_M
def test11_1():
    # Define the parameters for the photomultiplier setup
    n = 10              # number of cascades
    delta_U_n = 100     # change in voltage of the nth stage (in volts)
    U_n = 1000          # original voltage of the nth stage (in volts)
    M = 1000            # original gain of the photomultiplier

    # Calculate the change in gain
    delta_M = delta_gain(n, delta_U_n, U_n, M)

    # Print the result
    print(f"зависимость изменения коэф. усиления от изм. напр. питания(умножители): {delta_M:.2f}")
def test11():
    # Define the parameters for the photomultiplier setup
    n = 10              # number of cascades
    delta_U_n = 100     # change in voltage of the nth stage (in volts)
    U_n = 1000          # original voltage of the nth stage (in volts)
    M = 1000            # original gain of the photomultiplier

    # Calculate the change in gain
    delta_M = delta_gain(n, delta_U_n, U_n, M)

    # Print the result
    print(f"Change in gain: {delta_M:.2f}")

def max_integral_sensitivity(sigma, U_n, theta):
    """
    Computes the maximum integral sensitivity of a bolometer with an amplifier
    at R_H = R.
    
    Args:
        alpha (float): the temperature coefficient of resistance (in 1/K)
        sigma (float): the total thermal conductivity (in W/K)
        U_n (float): the supply voltage (in V)
    
    Returns:
        float: the maximum integral sensitivity (in V/W)
    """
    alpha = 1 / theta
    return (1/4) * (alpha / sigma) * U_n
def test12_1():
    # define the parameters
    sigma = 1.2
    U_n = 5.0
    theta = 0.02

    # calculate the maximum integral sensitivity
    S_v = max_integral_sensitivity(sigma, U_n, theta)

    # print the result
    print('Максимальная интегральная чувствительность = ', S_v)
def test12():
    # define the parameters
    sigma = 1.2
    U_n = 5.0
    theta = 0.02

    # calculate the maximum integral sensitivity
    S_v = max_integral_sensitivity(sigma, U_n, theta)

    # print the result
    print('maximum integral sensitivity = ', S_v)

# флуктуация потоков излучения фона и приемника
def noise_power(background_temp, receiver_temp, absorption_coeff, frequency_band, photodetector_area):
    k = 1.38e-23  # Boltzmann constant
    sigma = 5.67e-8  # Stefan-Boltzmann constant
    # Calculate the noise power due to background radiation
    background_power = 8 * k * sigma * photodetector_area * absorption_coeff * frequency_band * (background_temp ** 5)

    # Calculate the noise power due to receiver temperature
    receiver_power = 8 * k * sigma * photodetector_area * absorption_coeff * frequency_band * (receiver_temp ** 4)

    # Calculate the total noise power
    total_power = background_power + receiver_power

    return total_power
def test13_1():
    background_temp = 300  # K
    receiver_temp = 290  # K
    absorption_coeff = 0.8
    frequency_band = 1e6  # Hz
    photodetector_area = 1e-4  # m^2

    result = noise_power(background_temp, receiver_temp, absorption_coeff, frequency_band, photodetector_area)
    print('флуктуация потоков излучения фона и премника: ', result)
def test13():
    background_temp = 300  # K
    receiver_temp = 290  # K
    absorption_coeff = 0.8
    frequency_band = 1e6  # Hz
    photodetector_area = 1e-4  # m^2

    result = noise_power(background_temp, receiver_temp, absorption_coeff, frequency_band, photodetector_area)
    print(result)

# дисперсия напряжения рад. шума 
def calculate_radiation_noise_variance(S_intBack, ΔF_back, S_intReceiverRadiation, ΔF_receiverRadiation):
    V_radiation_squared = S_intBack**2 * ΔF_back**2 + S_intReceiverRadiation**2 * ΔF_receiverRadiation**2
    return V_radiation_squared
def test14():
    # Example usage
    S_intBack = 0.1  # integral sensitivity to background radiation
    ΔF_back = 0.5    # bandwidth of background radiation
    S_intReceiverRadiation = 0.2  # integral sensitivity to receiver radiation
    ΔF_receiverRadiation = 5    # bandwidth of receiver radiation

    V_radiation_squared = calculate_radiation_noise_variance(S_intBack, ΔF_back, S_intReceiverRadiation, ΔF_receiverRadiation)
    print("дисперсия напряжения рад. шума:", V_radiation_squared)


def calc_thermal_noise_voltage(T, R, bandwidth):
    k = 1.38e-23
    return 4 * k * T * R * bandwidth
def calc_thermal_noise_current(T, R, bandwidth):
    k = 1.38e-23
    G = 1 / R
    return 4 * k * T * G * bandwidth
def test15():
    # Calculate thermal noise voltage for a resistor with a resistance of 10 ohms, at a temperature of 300 K, over a bandwidth of 1 kHz
    voltage = calc_thermal_noise_voltage(300, 1000, 1000)
    print("Термический шум напряжения:", voltage, "V")

    # Calculate thermal noise current for a resistor with a resistance of 100 ohms, at a temperature of 273 K, over a bandwidth of 10 kHz
    current = calc_thermal_noise_current(273, 1000, 10000)
    print("Термический шум тока:", current, "A")


def current_shot_noise_power_density(I_0, delta_f):
    q = 1.6e-19  # Elementary charge in coulombs
    """
    Calculates the power spectral density of current shot noise.

    Parameters:
    -----------
    I_0 : float
        Average current in amperes.
    delta_f : float
        Bandwidth in hertz.

    Returns:
    --------
    P_i : float
        Power spectral density of current shot noise in watts per hertz.
    """
    return 2 * q * I_0 * delta_f
def rms_current_shot_noise(I_0, delta_f):
    """
    Calculates the RMS current of shot noise.

    Parameters:
    -----------
    I_0 : float
        Average current in amperes.
    delta_f : float
        Bandwidth in hertz.

    Returns:
    --------
    delta_I_shot : float
        RMS current of shot noise in amperes.
    """
    P_i = current_shot_noise_power_density(I_0, delta_f)
    return math.sqrt(P_i)
def test16():
    delta_f = 1e3  # frequency band in Hz
    I_0 = 1e-9  # average current in Amperes

    print(f"Дробовый шум тока: {rms_current_shot_noise(I_0, delta_f)} A")

def rms_current_shot_voltage(I_0, R_load, delta_f):
    e = 1.6e-19 # electron charge in coulombs

    return math.sqrt(2 * e * I_0 * R_load**2 * delta_f)
def test17():
    I_0 = 1e-6 # average current in amperes
    R_load = 200 # load resistance in ohms
    delta_f = 1e3 # frequency bandwidth in hertz

    V_shot = rms_current_shot_voltage(I_0, R_load, delta_f)
    print("Дробовое напряжение = ", V_shot)


def avalanche_shot_noise(I_DC, M, delta_f):
    """
    Calculates the shot noise current in an avalanche photodiode.

    Args:
        I_DC (float): DC current through the device in amperes.
        M (float): Gain of the avalanche process.
        delta_f (float): Frequency bandwidth of the measurement in hertz.

    Returns:
        float: Shot noise current in amperes.
    """
    q = 1.6e-19  # elementary charge in coulombs
    I_shot_squared = 2 * q * M ** 3 * I_DC * delta_f
    I_shot = I_shot_squared ** 0.5  # taking square root of variance to get standard deviation
    return I_shot
def test18():
    I_DC = 1e-6  # 10 microamps
    M = 100
    delta_f = 1e3  # 1 MHz
    I_shot = avalanche_shot_noise(I_DC, M, delta_f)
    print(f"Шум лавинного тока: {I_shot:.2e} A")

def voltage_gr_noise(V_supply, R_load, R_dark, R_T, tau_carrier, n, volume, f, delta_f):
    """
    Calculate the voltage dispersion of the generation-recombination noise for photoresistors.

    Args:
    V_supply (float): The supply voltage
    R_load (float): The load resistance
    R_dark (float): The dark resistance of the photoresistor
    tau_carrier (float): The lifetime of the carriers
    n (float): The concentration of carriers
    volume (float): The volume of the photo layer
    f (float): The frequency of modulation of the radiation flux
    delta_f (float): The frequency band

    Returns:
    float: The voltage dispersion of the generation-recombination noise for photoresistors
    """

    k = 1.38e-23  # Boltzmann's constant
    T = 300  # Room temperature in Kelvin
    q = 1.6e-19  # Elementary charge

    omega = 2 * math.pi * f

    # C = n * volume * q**2 / (2 * epsilon * k * T)

    term1 = (R_load**2 * R_dark**2) / (R_load + R_T)

    term2 = (tau_carrier * delta_f) / (n * volume)

    term3 = 1 / (1 + (2 * math.pi * f * tau_carrier)**2)

    return 4 * V_supply**2 * term1 * term2 * term3
def test19():
    V_gr = voltage_gr_noise(200, 1000, 100, 25, 1e-6, 1e18, 1e-9, 1000, 10000)
    print('Фоновый шум напряжения: ', V_gr)

def flicker_noise(I_0, B, A, f):
    e = 1.6e-19
    return 2 * e * I_0 * (1 + (B * I_0) / (A * f))
def test20():
    I_0 = 1e-9  # average value of current in the photoelement
    A = 1e-4   # area of the photocathode
    B = 1.5   # constant that depends on the photocathode
    f = 10e3  # frequency of modulation of the radiation flux

    print('Дисперсия фонового мерцания = ', flicker_noise(I_0, A, B, f))  # returns the value of the flicker noise dispersion

def voltage_dispersion_current_noise(R, I_0, B, delta_f, f_M):
    return (B * R**2 * I_0**2 * delta_f) / f_M
def test21():
    print('Дисперсия напряжения = ', voltage_dispersion_current_noise(R=100, I_0=1e-9, B=0.5, delta_f=1e6, f_M=10e6))

def total_dispersion():
    V_radiation_squared = calculate_radiation_noise_variance(S_intBack=0.1, ΔF_back=0.5, S_intReceiverRadiation=0.2, ΔF_receiverRadiation=5)
    voltage = calc_thermal_noise_voltage(300, 10, 1000)
    V_shot = rms_current_shot_voltage(I_0=1e-6, R_load=100, delta_f=1e6)
    V_gr = voltage_gr_noise(5, 1000, 100, 25, 1e-6, 1e18, 1e-9, 1000, 10000)
    V_cn = voltage_dispersion_current_noise(R=100, I_0=1e-9, B=0.5, delta_f=1e6, f_M=10e6)

    return V_radiation_squared**2 + voltage**2 + V_shot**2 + V_gr**2 + V_cn**2
def test22():
    print('Итоговая дисперсия = ', total_dispersion())

def rmsAmplyfier(U_out, K):
    return math.sqrt(U_out) / K
def test23():
    print('Среднеквадратическое усиления = ', rmsAmplyfier(U_out=5, K=3))

def ThreshholdFlow(S_U):
    return math.sqrt(total_dispersion()) / S_U
def test24():
    print('Пороговый поток = ', ThreshholdFlow(75))


# Define example functions for eps, tau, and M
def eps_func(lambda_val, Theta):
    return lambda_val + Theta
def tau_func(lambda_val):
    return lambda_val
def M_func(lambda_val, Theta):
    return lambda_val - Theta
def radiation_flux(q, D_in, f, eps_func, tau_func, M_func, lambda1, lambda2, Theta):
    integrand = lambda lambda_var: eps_func(lambda_var, Theta) * tau_func(lambda_var) * M_func(lambda_var, Theta)
    F = (q / 4) * (D_in / f)**2 * quad(integrand, lambda1, lambda2)[0]
    return F
def test25():
    # Define values for the parameters
    q = pi * (0.1 / 2)**2  # Area of target aperture
    D_in = 0.1  # Diameter of entrance pupil
    f = 1.0  # Focal length of the lens
    # Define integration limits
    lambda1 = 40e-2
    lambda2 = 70e-2

    # Calculate radiation flux
    F = radiation_flux(q, D_in, f, eps_func, tau_func, M_func, lambda1, lambda2, Theta=0.5)

    # Print result
    print('СВетовой поток = ', F) # Output: 0.00023713873066238818
    
def current_photodetector_flux(q, D_in, f, S_func, eps_func, tau_func, M_func, lambda1, lambda2, Theta):
    integrand = lambda x: S_func(x) * eps_func(x, Theta) * tau_func(x) * M_func(x, Theta)
    F = (q/4) * (D_in/f)**2 * quad(integrand, lambda1, lambda2)[0]
    return F
def test26():
    # Example usage
    q = 1.5 # cm^2
    D_in = 2 # cm
    f = 10 # cm
    lambda1 = 400 # nm
    lambda2 = 700 # nm
    Theta = 0 # degrees
    S_func = lambda x: 0.5 # some arbitrary function for demonstration purposes
    eps_func = lambda x, Theta: 0.8 # some arbitrary function for demonstration purposes
    tau_func = lambda x: 0.9 # some arbitrary function for demonstration purposes
    M_func = lambda x, Theta: 1.0 # some arbitrary function for demonstration purposes

    flux = current_photodetector_flux(q, D_in, f, S_func, eps_func, tau_func, M_func, lambda1, lambda2, Theta)
    print('Ток фотоприемника = ', flux)

def s_lambda(x):
    return np.exp(-((x - 550) / 50)**2)
def absolute_sensitivity(S, s_lambda, K_u):
    """
    Computes the absolute sensitivity given the integral sensitivity S, the relative spectral sensitivity s_lambda,
    and the utilization factor for the reference emitter K_u.
    
    Args:
    - S (float): the integral sensitivity of the photodetector
    - s_lambda (function): a function that takes a wavelength as input (in nm) and returns the relative spectral 
                           sensitivity of the photodetector at that wavelength
    - K_u (float): the utilization factor for the reference emitter
    
    Returns:
    - S_lambda (function): a function that takes a wavelength as input (in nm) and returns the absolute spectral 
                           sensitivity of the photodetector at that wavelength
    
    """
    S_lambda_max = S * K_u
    S_lambda = lambda x: S_lambda_max * s_lambda(x)
    
    return S_lambda
def test27():
    S_lambda = absolute_sensitivity(1, s_lambda, 0.8)

    # Evaluate the absolute sensitivity at a few wavelengths
    print('Сравнение абсолютной чувствительности при 500, 550, 600 длинах волн\n')
    print(S_lambda(500))  # Output: 0.15665674855704327
    print(S_lambda(550))  # Output: 0.8000000000000002
    print(S_lambda(600))  # Output: 0.1953668148131644
    # Define the integral sensitivity S, and the utilization factor K_u

# Define the necessary functions and constants
def S_lambda(s_lambda_max, lambda_min, lambda_max, lambda_var):
    # Implementation of S_lambda function
    return 0.8
def eps(lambda_var, theta):
    # Implementation of eps function
    return lambda_var + theta
def tau(lambda_var):
    # Implementation of tau function
    return lambda_var
def M_lambda_0(lambda_var, theta):
    # Implementation of M_lambda_0 function
    return lambda_var - theta
def A(q, D_in, f, S_lambda_max, R_l, R_i, K_g):
    return (q / 4) * (D_in / f)**2 * ((S_lambda_max * R_l * R_i * K_g) / (R_l + R_i))
def photocurrent_voltage(lambda_min, lambda_max, theta, s_lambda_max, q, D_in, f, R_l, R_i, K_g):
    # Calculate the necessary integrand
    integrand = lambda lambda_var: S_lambda(s_lambda_max, lambda_min, lambda_max, lambda_var) * eps(lambda_var, theta) * tau(lambda_var) * M_lambda_0(lambda_var, theta)

    # Calculate A
    A_val = A(q, D_in, f, s_lambda_max, R_l, R_i, K_g)

    # Calculate the photocurrent voltage
    return A_val * quad(integrand, lambda_min, lambda_max)[0]
def test28():
    lambda_min = 400 # nm
    lambda_max = 700 # nm
    theta = 45 # degrees
    s_lambda_max = 1.0 # arbitrary units
    q = 1.6e-19 # C
    D_in = 5e-3 # m
    f = 10e-3 # m
    R_l = 1e3 # ohm
    R_i = 1e6 # ohm
    K_g = 1e6 # V/A

    voltage = photocurrent_voltage(lambda_min, lambda_max, theta, s_lambda_max, q, D_in, f, R_l, R_i, K_g)

    print('Напряжение фотока = ', voltage)

# DON'T WORK!!!!!!!!!!!
def photocurrent_voltage_rectangle(lambda_min, lambda_max, theta, s_lambda_max, q, D_in, f, R_l, R_i, K_g, n):
    # Calculate the necessary integrand
    integrand = lambda lambda_var: S_lambda(s_lambda_max, lambda_min, lambda_max, lambda_var) * eps(lambda_var, theta) * tau(lambda_var) * M_lambda_0(lambda_var, theta)

    # Calculate A
    A_val = A(q, D_in, f, s_lambda_max, R_l, R_i, K_g)

    # Calculate the step size
    delta_lambda = (lambda_max - lambda_min) / n

    # Initialize the sum
    integral_sum = 0

    # Iterate over the intervals
    for i in range(n):
        # Calculate the midpoint of the interval
        lambda_mid = lambda_min + (i + 0.5) * delta_lambda

        # Add the integrand multiplied by the step size to the sum
        integral_sum += integrand(lambda_mid) * delta_lambda

    # Calculate the photocurrent voltage
    return A_val * s_lambda_max * integral_sum
def test29():
    # Define the required parameters
    lambda_min = 300 # nm
    lambda_max = 800 # nm
    theta = 0 # degrees
    s_lambda_max = 0.8
    q = 1.6e-19 # C
    D_in = 1e-3 # m
    f = 10 # Hz
    R_l = 1e3 # ohm
    R_i = 1e3 # ohm
    K_g = 1e3

    # Calculate the photocurrent voltage
    result = photocurrent_voltage_rectangle(lambda_min, lambda_max, theta, s_lambda_max, q, D_in, f, R_l, R_i, K_g, 100)

    print(result)


def numerical_rectangle(A, S_max, s_func, eps_func, tau_func, M_func, lambda1, lambda2, Theta, delta_lambda):
    # Find the number of sampling points
    N = int((lambda2 - lambda1) / delta_lambda)
    
    # Calculate the sum
    total_sum = 0
    for i in range(N):
        # Calculate the current coordinate
        lambda_i = lambda1 + i * delta_lambda
        # Calculate the discretized function
        M_0_lambda_i_Theta = M_func(lambda_i, Theta)
        # Calculate the integrand at the current point
        integrand = s_func(lambda_i) * eps_func(lambda_i, Theta) * tau_func(lambda_i) * M_0_lambda_i_Theta
        # Add the integrand to the total sum
        total_sum += integrand
    
    # Multiply the sum by the appropriate constants
    result = A * S_max * total_sum * delta_lambda
    
    return result
def test30():
    # Define the parameters
    # Define the required parameters
    s_lambda_max = 0.8
    q = 1.6e-19 # C
    D_in = 1e-3 # m
    f = 10 # Hz
    R_l = 1e3 # ohm
    R_i = 1e3 # ohm
    K_g = 1e3

    A_val = A(q, D_in, f, s_lambda_max, R_l, R_i, K_g)
    S_max = 0.8
    s_func = lambda x: x
    eps_func = lambda x, y: x + y
    tau_func = lambda x: x
    M_func = lambda x, y: x + y
    lambda1 = 300
    lambda2 = 800
    Theta = 0.0
    delta_lambda = 0.1

    # Approximate the value of U_f
    result = numerical_rectangle(A_val, S_max, s_func, eps_func, tau_func, M_func, lambda1, lambda2, Theta, delta_lambda)

    # Print the result
    print('Нормализация фототока прямоугольниками = ', result)

def D_sum_in(D, K):
    """
    Calculates the output variance of a system with component variances D and scaling factors K,
    assuming no correlations between the component errors.
    
    Args:
    - D (list or tuple of floats): Component variances
    - K (list or tuple of floats): Scaling factors
    
    Returns:
    - float: Output variance
    """
    n = len(D)
    D_sum_in = 0.0
    for i in range(n):
        D_sum_in += (K[i]**2) * (D[i] * (K[i+1]**2) if i < n-1 else D[i])
    return D_sum_in
def test31():
    D = [0.1, 0.2, 0.3]
    K = [1.5, 2.0, 1.8]

    D = D_sum_in(D, K)

    print('Выходная дисперсия системы = ', D) # prints the total output variance of the system


def fullDispersion(D, K):
    K_full = 1
    for i in K:
        K_full *= i
    return D_sum_in(D, K) * K_full**-2
def test32():
    D = [0.1, 0.2, 0.3]
    K = [1.5, 2.0, 1.8]

    Di = fullDispersion(D, K)

    print('Итоговая дисперсия = ', Di) # prints the total output variance of the system

def unidirectional_scanning_error(velocity_product, scanning_speed):
    return velocity_product / (scanning_speed - velocity_product)
def test33():
    velocity_product = 10  # meters per second
    scanning_speed = 30  # meters per second
    error = unidirectional_scanning_error(velocity_product, scanning_speed)
    print('Погрешность однонаправленного сканирования', error)  # Output: 0.5

def alternating_scanning_error(velocity_product, scanning_speed):
    return velocity_product**2 / (scanning_speed**2 - velocity_product**2)
def test34():
    velocity_product = 10  # meters per second
    scanning_speed = 30  # meters per second
    error = alternating_scanning_error(velocity_product, scanning_speed)
    print('Погрешность двунаправленного сканирования = ', error)  # Output: 0.5


def unidirectional_scanning_speed(vs, sigma1):
    return vs * (1 + 1 / sigma1)
def bidirectional_scanning_speed(vs, sigma2):
    return vs * math.sqrt(1 + 1 / sigma2)
def test35():
    print('скорость однонаправленного сканирования = ', unidirectional_scanning_speed(vs=10, sigma1=unidirectional_scanning_error(10, 30)))

    print('скорость двунаправленного сканирования = ', bidirectional_scanning_speed(vs=10, sigma2=alternating_scanning_error(10, 30)))


def calculate_information_per_scan(S_sc, delta, K_m):
    N = math.log2((S_sc * K_m) / delta)
    return N
def calculate_information_per_scan_secod(S_sc_prime, delta, K_m):
    N = math.log2(S_sc_prime / delta)
    return N
def test36():
    S_sc = 10 # mm
    S_sc_prime = 5 # mm
    delta = 0.1 # mm
    K_m = 20

    N = calculate_information_per_scan(S_sc_prime, delta, K_m)
    print(f"Количество информации за такт однонаправленного сканирования: {N:.2f} бит")

    N = calculate_information_per_scan_secod(S_sc, delta, K_m)
    print(f"Количество информации за такт двунаправленного сканирования: {N:.2f} бит")


def calc_info_per_scan(S_sc_prime, K_m, delta_s):
    Delta_s_prime = delta_s / K_m
    N_prime = math.log2(S_sc_prime / Delta_s_prime)
    return N_prime
def test37():
    N_prime = calc_info_per_scan(S_sc_prime=5, K_m=20, delta_s=0.1)
    print('Количество информации за такт скана = ', N_prime)


def calculate_spectral_efficiency(s, Φ_eλ, Φ_eλ_max):
    """
    Calculates the spectral efficiency coefficient of the system based on the
    spectral distribution of the input radiation flux, the relative spectral
    characteristic of the transmission or transformation of the radiation flux,
    and the maximum spectral density of the input radiation flux.

    Args:
    s (array-like): Array of values representing the relative spectral
                    characteristic of the transmission or transformation of
                    the radiation flux.
    Φ_eλ (array-like): Array of values representing the spectral distribution
                       of the spectral density of the input radiation flux.
    Φ_eλ_max (float): Maximum spectral density of the input radiation flux.

    Returns:
    float: Spectral efficiency coefficient of the system.
    """
    numerator = np.trapz(s * Φ_eλ, dx=1)
    denominator = np.trapz(Φ_eλ / Φ_eλ_max, dx=1)
    return numerator / denominator
def test38():
    # Example input arrays and max spectral density
    s = np.array([0.8, 0.9, 0.95, 1.0, 0.95, 0.9, 0.8])
    Φ_eλ = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6])
    Φ_eλ_max = 1.0

    # Calculate spectral efficiency coefficient
    χ = calculate_spectral_efficiency(s, Φ_eλ, Φ_eλ_max)

    # Print result
    print(f"коэффициент спектральной эффективности системы: {χ}")
#OR ALTERNATIVE
# Define the s_lambda and Phi_e_lambda functions
def s_lambda(lambda_):
    # Define the relative spectral characteristic of the transmission or transformation of the radiation flux
    # For example, let's say s_lambda is a Gaussian function centered at 500 nm with a standard deviation of 10 nm
    center_lambda = 500  # nm
    sigma_lambda = 10  # nm
    return np.exp(-(lambda_ - center_lambda)**2 / (2 * sigma_lambda**2))
def Phi_e_lambda(lambda_):
    # Define the spectral distribution of the spectral density of the input radiation flux
    # For example, let's say Phi_e_lambda is a constant function equal to 1 between 400 nm and 600 nm, and 0 elsewhere
    if 400 <= lambda_ <= 600:
        return 1
    else:
        return 0
def spectral_efficiency_coefficient(s_lambda, Phi_e_lambda, max_lambda):
    numerator_integrand = lambda lambda_: s_lambda(lambda_) * Phi_e_lambda(lambda_)
    denominator_integrand = lambda lambda_: Phi_e_lambda(lambda_)
    
    numerator = quad(numerator_integrand, 0, max_lambda)[0]
    denominator = quad(denominator_integrand, 0, max_lambda)[0]
    
    return numerator / denominator
def test39():
    # Define the maximum value of lambda
    max_lambda = 1000  # nm

    # Calculate the spectral efficiency coefficient
    chi = spectral_efficiency_coefficient(s_lambda, Phi_e_lambda, max_lambda)
    print('То же но альтернативным методом = ', chi)

# КРАЙНЕ ВАЖНО КРАЙНЕ ВАЖНО КРАЙНЕ ВАЖНО!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def planck_law(wavelength, T):
    h = 6.625e-34  # Planck constant [J s]
    c = 2.998e+8   # speed of light [m/s]
    k = 1.381e-23  # Boltzmann constant [J/K]

    c1 = 2 * np.pi * h * np.power(c, 2)
    c2 = h * c / k
    expon = c2 / (wavelength * T)
    res = c1 * np.power(wavelength, -5.0) * np.power((np.exp2(expon) - 1.0), -1.0)
    return res
def test40():
    wavelength = np.linspace(1e-9, 3e-6, 1000)  # wavelength range [m]
    T = 300                                    # temperature [K]
    spectrum = planck_law(wavelength, T)       # spectral energy density [W/m^3]
    print('Закон Планка, тестирование = ', spectrum)
def planck_spectral_density(wavelength, temperature):
    h = 6.625 * 10**-34
    c = 2.98 * 10**10
    k = 1.3805 * 10**-23
    c1 = 3.74 * 10**-12
    c2 = 1.438 * 10**4
    expon = c2 / (wavelength * temperature)
    if expon > 700:
        return 0.0
    else:
        return c1 * np.power(wavelength, -5.0) * np.power((np.exp(expon) - 1.0), -1.0)
def test40_2():
    wavelength = 500e-9  # 500 nm
    temperature = 300  # Kelvin
    spectral_density = planck_spectral_density(wavelength, temperature)
    print('То же но другим методом = ', spectral_density)


def lamberts_law(B, dS, alpha):
    return B * dS * np.cos(alpha)
def test41():
    brightness = 100 # W/m^2/sr
    radiating_area = 0.5 # m^2
    angle = 45 # degrees

    intensity = lamberts_law(brightness, radiating_area, angle)
    print("Интенсивность света под углом {0} градусов составляет {1:.2f} Вт/м^2".format(angle, intensity))



def transformed_planck_law(wavelength, temperature):
    # Constants
    h = 6.626e-34  # Planck constant
    c = 2.998e8  # Speed of light
    k = 1.381e-23  # Boltzmann constant

    c1 = 2 * np.pi * h * c**2
    c2 = h * c / k

    x = wavelength * temperature / c2
    numerator = c2**5 * x**(-5)
    denominator = c1 * (np.exp(x) - 1.0)
    return numerator / denominator
def spectral_radiance(wavelength, temperature, solid_angle):
    return transformed_planck_law(wavelength, temperature) * np.pi * solid_angle
def test42():
    wavelength = 600e-9  # 500 nm
    temperature = 300  # Kelvin
    solid_angle = 45 * np.pi * (1 - np.cos(np.deg2rad(1)))  # steradians
    radiance = spectral_radiance(wavelength, temperature, solid_angle)
    print(f"(ПЛОХО РАБОТАЕТ)Спектральное излучение: {radiance:.4f} Вт/(м^2 ср·нм)")

def planck_transformed(x):
    h = 6.626e-34
    c = 2.998e8
    k = 1.381e-23
    c1 = 2 * np.pi * h * c**2
    c2 = h * c / k
    r = c2**5 / (c1 * x**5)
    y = r / (np.exp(1/x) - 1)
    return y / np.max(y)  # normalization
def test43():
    # Plot the transformed Planck's law for a range of x values
    x_vals = np.linspace(1e-10, 5, 10)  # adjust the range as needed
    y_vals = planck_transformed(x_vals)

    plt.plot(x_vals, y_vals)
    plt.xlabel('x(длины)')
    plt.ylabel('y(нормализация спектральной хар-ки)')
    plt.title('Спектральная характеристика по Планку')
    plt.show()


def transmittance_lambda(lambda_):
    # define your own transmittance function here
    return 0.5 # example transmittance function that always returns 0.5
def spectral_efficiency(transmittance_lambda, max_lambda, theta):
    c1 = 3.74183e-12 # in W.m^2
    c2 = 1.4388 # in cm.K
    x_max = theta / c2
    y_integrand = lambda lambda_: 142.32 * np.power((lambda_ / max_lambda), -5.0) * np.power(np.exp(4.9651 / (lambda_ * theta / c2)) - 1.0, -1.0)
    numerator_integrand = lambda lambda_: transmittance_lambda(lambda_) * y_integrand(lambda_)
    denominator_integrand = lambda lambda_: y_integrand(lambda_)
    
    numerator = quad(numerator_integrand, 0, np.inf)[0]
    denominator = quad(denominator_integrand, 0, np.inf)[0]
    
    return numerator / denominator
def test44():
    max_lambda = 1000 # example maximum wavelength value
    theta = 5000 # example temperature value

    print(spectral_efficiency(transmittance_lambda, max_lambda, theta)) # calculate spectral efficiency with example inputs



def s(wavelength):
    # define the relative spectral sensitivity function of the photodetector
    # this is just an example, replace with the appropriate function for your detector
    return 0.8

def y_1(wavelength, temperature):
    # define the dimensionless value of the spectral density function of the energy luminosity at the output of the optical element
    # this is just an example, replace with the appropriate function for your system
    return 142.32 * ((wavelength * temperature / 1.4388e-2)**-5.0) * ((np.exp(1.4388e-2 / (wavelength * temperature)) - 1.0)**-1.0)

def photodetector_spectral_efficiency(s, y_1, wavelength_range):
    """
    Calculates the spectral efficiency of a photodetector.

    Parameters:
    s (callable): A function that takes a wavelength as input and returns the relative spectral sensitivity of the photodetector.
    y_1 (callable): A function that takes a wavelength and temperature as input and returns the dimensionless value of the spectral density function of the energy luminosity at the output of the optical element.
    wavelength_range (tuple): A tuple containing the lower and upper limits of the wavelength range over which to integrate.

    Returns:
    The spectral efficiency of the photodetector.
    """
    wavelength_lower, wavelength_upper = wavelength_range
    integral_num = lambda x, temperature: s(x) * y_1(x, temperature)  # Numerator integrand
    integral_denom = y_1  # Denominator integrand
    numerator, _ = quad(integral_num, wavelength_lower, wavelength_upper, args=(125,))
    denominator, _ = quad(integral_denom, wavelength_lower, wavelength_upper, args=(125,))
    return numerator / denominator


def test45():
    wavelength_range = (1e-9, 2e-6) # define the wavelength range over which to integrate
    temperature = 125
    spectral_efficiency = photodetector_spectral_efficiency(s, y_1, wavelength_range)
    print(spectral_efficiency)

#chat?hist=y94zr6wPtHiy-VmtKsA3a7qLFF5Ltpo-p0LTi5RmRtk