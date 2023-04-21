import numpy as np
import math
from scipy.integrate import quad
from typing import Callable
from scipy.constants import h, c, k
from math import pi
import matplotlib.pyplot as plt
from tests import *

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

sensitivity = photoresistor_sensitivity(C, U, a, F, b)
print("Средняя чувствительность:", sensitivity, "А/Вт/м^2")

differential_sensitivity = photoresistor_differential_sensitivity(C, U, a, F, b)
print("Дифференцированная чувствительность:", differential_sensitivity, "А/Вт/м^2")

V = 5
T = 700  # Temperature in Kelvin
I_sat = 1e-3  # Saturation current in Amperes
A = 4   # Coefficient for the material of the photodetector

dark_current = photoresistor_dark_current(V, T, I_sat, A)*-1
print("Темновой ток:", dark_current, "A")

trans = photodiode_transformations(C, U, a, F, b, T, I_sat, A)
print('преобразование фотодиода: ', trans)

photVolt = photoresistor_voltage(F, A, T, I_sat, C, U, a, b)
print('зависимость напряжения холостого хода на фотодиоде: ', photVolt)

def calculate_Un(wavelength):
    h = 6.62607015e-34   # Planck constant in J*s
    c = 299792458        # Speed of light in m/s
    q = 1.602176634e-19  # Elementary charge in C

    return (h * c / q) / wavelength
Sl = 0.9  # Amperes per Watt
Un = calculate_Un(500e-9)
U_max = calculate_max_product(sensitivity, Sl, Un, I_sat) / 10
print('Максимальное напряжение в фотодиодном режиме: ', U_max)


R_0 = 1000  # Ohms
I_f = 0.01  # Amperes
result = calc_S_Umax(sensitivity, R_0, Sl, I_f, dark_current)*-1
print("В фотогальваническом режиме = {:.4f} Вольт".format(result))

n = 10              # number of cascades
delta_U_n = 10     # change in voltage of the nth stage (in volts)
U_n = 100          # original voltage of the nth stage (in volts)
M = 100            # original gain of the photomultiplier
delta_M = delta_gain(n, delta_U_n, U_n, M)
print(f"зависимость изменения коэф. усиления от изм. напр. питания(умножители): {delta_M:.2f}")

sigma = 0.5
S_v = max_integral_sensitivity(sigma, Un, T)
print('Максимальная интегральная чувствительность = ', S_v)

receiver_temp = 290  # K
absorption_coeff = 0.8
frequency_band = 1e4  # Hz
photodetector_area = 1e-4  # m^2
result = noise_power(T, receiver_temp, absorption_coeff, frequency_band, photodetector_area)
print('флуктуация потоков излучения фона и премника: ', result)


ΔF_back = 0.5    # bandwidth of background radiation
S_intReceiverRadiation = 0.02  # integral sensitivity to receiver radiation
ΔF_receiverRadiation = 5    # bandwidth of receiver radiation
V_radiation_squared = calculate_radiation_noise_variance(S_v, ΔF_back, S_intReceiverRadiation, ΔF_receiverRadiation)
print("дисперсия напряжения рад. шума:", V_radiation_squared)


voltage2 = calc_thermal_noise_voltage(T, R_0, 1000)
print("Термический шум напряжения:", voltage2, "V")
current = calc_thermal_noise_current(T, R_0, 10000)
print("Термический шум тока:", current, "A")


print(f"Дробовый шум тока: {rms_current_shot_noise(current, frequency_band)} A")


R_load = 200 # load resistance in ohms
V_shot = rms_current_shot_voltage(current, R_load, frequency_band)
print("Дробовое напряжение = ", V_shot)

M = 100
I_shot = avalanche_shot_noise(dark_current, M, frequency_band)
print(f"Шум лавинного тока: {I_shot:.2e} A")


V_gr = voltage_gr_noise(5, R_load, 100, 25, 1e-6, 1e18, 1e-9, 1000, frequency_band)
print('Фоновый шум напряжения: ', V_gr)

I_0 = 0.001  # average value of current in the photoelement
A = 1e-4   # area of the photocathode
B = 1.5   # constant that depends on the photocathode
f = 10e3  # frequency of modulation of the radiation flux
print('Дисперсия фонового мерцания = ', flicker_noise(I_0, A, B, f))  # returns the value of the flicker noise dispersion


print('Дисперсия напряжения = ', voltage_dispersion_current_noise(R_load, I_0, B, frequency_band, f_M=10e6))

print('Итоговая дисперсия = ', 0.010000784289870465 + 3.8640000000000006e-14 + 7.032723512267491e-15 + 0.0017777075966948967 + 6e-05)

print('Среднеквадратическое усиления = ', rmsAmplyfier(U_out=5, K=3))

print('Пороговый поток = ', ThreshholdFlow(75))

q = pi * (0.1 / 2)**2  # Area of target aperture
D_in = 0.1  # Diameter of entrance pupil
f = 1.0  # Focal length of the lens
lambda1 = 40e-2
lambda2 = 70e-2
F = radiation_flux(q, D_in, f, eps_func, tau_func, M_func, lambda1, lambda2, Theta=0.5)
print('Световой поток = ', F) # Output: 0.00023713873066238818