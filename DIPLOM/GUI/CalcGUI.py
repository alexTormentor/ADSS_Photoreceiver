from tkinter import *
from tests import *
import ADSS_FinalVer

SMTH = ADSS_FinalVer.CalcSelector2

class PhoRes(SMTH):
    def __init__(self, master):
        self.C = DoubleVar()
        self.U = DoubleVar()
        self.W = DoubleVar()
        self.Time = DoubleVar()
        self.master = master

        self.master.geometry('600x300+100+200')
        self.master.title("ещё одно окошко =)")
        self.label1 = Label(self.master, text="Введите параметры", fg='red').grid(row=0, column=2)
        self.Clabel = Label(self.master, text='Введите константу').grid(row=1, column=0)
        self.Ulabel = Label(self.master, text='Введите напряжение').grid(row=2, column=0)
        self.Wlabel = Label(self.master, text='Введите энергию').grid(row=3, column=0)
        self.Timelabel = Label(self.master, text='Введите время').grid(row=4, column=0)
        self.CEntry = Entry(self.master, textvar=self.C).grid(row=1, column=2)
        self.UEntry = Entry(self.master, textvar=self.U).grid(row=2, column=2)
        self.WEntry = Entry(self.master, textvar=self.W).grid(row=3, column=2)
        self.TimeEntry = Entry(self.master, textvar=self.Time).grid(row=4, column=2)
        self.button1 = Button(self.master, text='рассчитать фототок', fg='blue',
                              command=self.photocur).grid(row=9, column=2)
        self.button2 = Button(self.master, text='ср. чувств-ь', fg='blue', command=self.avSens).grid(row=9, column=3)
        self.button3 = Button(self.master, text='дифф. чувств-ь', fg='blue',
                              command=self.difSens).grid(row=10, column=3)
        self.button0 = Button(self.master, text='Назад', fg='blue', command=self.quit).grid(row=5, column=3)



    def difSens(self):
        global state
        C = self.C.get()
        U = self.U.get()
        W = self.W.get()
        Time = self.Time.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        differential_sensitivity = photoresistor_differential_sensitivity(C, U, a, F, b)
        print(differential_sensitivity)
        state = "Дифференцированная чувствительность =  " + str(differential_sensitivity)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=2)
        return state

    def avSens(self):
        global state
        C = self.C.get()
        U = self.U.get()
        W = self.W.get()
        Time = self.Time.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        sensitivity = photoresistor_sensitivity(C, U, a, F, b)
        print(sensitivity)
        state = "Средняя чувствительность =  " + str(sensitivity)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=2)
        return state

    def photocur(self):
        global state
        C = self.C.get()
        U = self.U.get()
        W = self.W.get()
        Time = self.Time.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        current2 = photoresistor_current(C, U, a, F, b)
        print(current2)
        state = "Ваш фототок =  " + str(current2)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=2)
        return state

    def quit(self):
        self.master.destroy()



class PhoDio(SMTH):
    def __init__(self, master):

        self.V = DoubleVar()
        self.T = DoubleVar()
        self.ISat = DoubleVar()
        self.A = DoubleVar()
        self.C = DoubleVar()
        self.E = DoubleVar()
        self.Time = DoubleVar()
        self.Sl = DoubleVar()
        self.Vn = DoubleVar()
        self.sens = DoubleVar()
        self.R = DoubleVar()
        self.darkCur = DoubleVar()
        self.I = DoubleVar()
        self.master = master

        self.master.geometry('1080x720+100+200')
        self.master.title("ещё одно окошко =)")
        self.label1 = Label(self.master, text="Введите параметры", fg='red').grid(row=0, column=2)
        self.Vlabel = Label(self.master, text='напряжение').grid(row=1, column=0)
        self.Tlabel = Label(self.master, text='температура').grid(row=2, column=0)
        self.ISatlabel = Label(self.master, text='ток насыщения').grid(row=3, column=0)
        self.Alabel = Label(self.master, text='тип диода').grid(row=4, column=0)
        self.Clabel = Label(self.master, text='коэффициент').grid(row=5, column=0)
        self.Elabel = Label(self.master, text='энергия').grid(row=6, column=0)
        self.Timelabel = Label(self.master, text='время').grid(row=7, column=0)
        self.Sllabel = Label(self.master, text='А/Вт').grid(row=8, column=0)
        self.Vnlabel = Label(self.master, text='длина волны').grid(row=9, column=0)
        self.Rlabel = Label(self.master, text='Сопротивление').grid(row=10, column=0)
        self.darkCurlabel = Label(self.master, text='темн. ток').grid(row=12, column=0)
        self.Ilabel = Label(self.master, text='ток').grid(row=13, column=0)
        self.VEntry = Entry(self.master, textvar=self.V).grid(row=1, column=2)
        self.TEntry = Entry(self.master, textvar=self.T).grid(row=2, column=2)
        self.ISatEntry = Entry(self.master, textvar=self.ISat).grid(row=3, column=2)
        self.AEntry = Entry(self.master, textvar=self.A).grid(row=4, column=2)
        self.CEntry = Entry(self.master, textvar=self.C).grid(row=5, column=2)
        self.EEntry = Entry(self.master, textvar=self.E).grid(row=6, column=2)
        self.TimeEntry = Entry(self.master, textvar=self.Time).grid(row=7, column=2)
        self.SlEntry = Entry(self.master, textvar=self.Sl).grid(row=8, column=2)
        self.VnEntry = Entry(self.master, textvar=self.Vn).grid(row=9, column=2)
        self.SensEntry = Entry(self.master, textvar=self.sens).grid(row=10, column=2)
        self.REntry = Entry(self.master, textvar=self.R).grid(row=11, column=2)
        self.darkCurEntry = Entry(self.master, textvar=self.darkCur).grid(row=12, column=2)
        self.IEntry = Entry(self.master, textvar=self.I).grid(row=13, column=2)

        self.button1 = Button(self.master, text='темн. ток', fg='blue', command=self.darkCurrent).grid(row=9, column=2)
        self.button2 = Button(self.master, text='преобр. диода', fg='blue', command=self.tranCur).grid(row=9, column=3)
        self.button3 = Button(self.master, text='зав. хода', fg='blue', command=self.turnDiode).grid(row=10, column=3)
        self.button4 = Button(self.master, text='режим диода', fg='blue', command=self.Vdiode).grid(row=10, column=4)
        self.button4 = Button(self.master, text='режим гальван', fg='blue', command=self.Vhalvan).grid(row=11, column=4)
        self.button0 = Button(self.master, text='Назад', fg='blue', command=self.quit).grid(row=5, column=3)



    def Vhalvan(self):
        global state
        R_0 = self.R.get()
        I_f = self.I.get()
        Sl = self.Sl.get()
        sens = self.sens.get()
        dark_cur = self.darkCur.get()


        result = calc_S_Umax(sens, R_0, Sl, I_f, dark_cur) * -1
        print(result)
        state = "В фотогальваническом режиме = {:.4f} Вольт " + str(result)
        self.labelresult = Label(self.master, text=state).grid(row=12, column=2)
        return state

    def Vdiode(self):
        global state
        I_sat = self.ISat.get()
        sensitivity = self.sens.get()
        Sl = self.Sl.get()

        def calculate_Un(wavelength):
            h = 6.62607015e-34  # Planck constant in J*s
            c = 299792458  # Speed of light in m/s
            q = 1.602176634e-19  # Elementary charge in C

            return (h * c / q) / wavelength


        Vn = calculate_Un(500e-9)


        U_max = calculate_max_product(1, 1, 1, 1) / 10
        print(U_max)
        state = "Максимальное напряжение в фотодиодном режиме =  " + str(U_max)
        self.labelresult = Label(self.master, text=state).grid(row=9, column=2)
        return state

    def turnDiode(self):
        global state
        C = self.C.get()
        E = self.E.get()
        Time = self.Time.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=E, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        V = self.V.get()
        T = self.T.get()
        I_sat = self.ISat.get()
        A = self.A.get()

        photVolt = photoresistor_voltage(F, A, T, I_sat, C, V, a, b)
        print(photVolt)
        state = "зависимость напряжения холостого хода на фотодиоде =  " + str(photVolt)
        self.labelresult = Label(self.master, text=state).grid(row=9, column=2)
        return state

    def tranCur(self):
        global state
        C = self.C.get()
        E = self.E.get()
        Time = self.Time.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=E, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        V = self.V.get()
        T = self.T.get()
        I_sat = self.ISat.get()
        A = self.A.get()

        trans = photodiode_transformations(C, V, a, F, b, T, I_sat, A)
        print(trans)
        state = "Преобразование диода =  " + str(trans)
        self.labelresult = Label(self.master, text=state).grid(row=10, column=2)
        return state


    def darkCurrent(self):
        global state
        V = self.V.get()
        T = self.T.get()
        I_sat = self.ISat.get()
        A = self.A.get()
        dark_current = photoresistor_dark_current(V, T, I_sat, A) * -1
        print(dark_current)
        state = "Темновой ток =  " + str(dark_current)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=2)
        return state


    def quit(self):
        self.master.destroy()


class PhoMul(SMTH):
    def __init__(self, master):
        self.n = DoubleVar()
        self.deltaUn = DoubleVar()
        self.Un = DoubleVar()
        self.M = DoubleVar()
        self.sigma = DoubleVar()
        self.wave = DoubleVar()
        self.T = DoubleVar()
        self.receiver_temp = DoubleVar()
        self.absorption_coeff = DoubleVar()
        self.frequency_band = DoubleVar()
        self.photodetector_area = DoubleVar()
        self.ΔF_back = DoubleVar()
        self.S_intReceiverRadiation = DoubleVar()
        self.ΔF_receiverRadiation = DoubleVar()
        self.S_v = DoubleVar()
        self.R0 = DoubleVar()
        self.I = DoubleVar()
        self.RLoad = DoubleVar()
        self.darkCur = DoubleVar()
        self.V_supply = DoubleVar()
        self.R_dark = DoubleVar()
        self.tau_carrier = DoubleVar()
        self.nCar = DoubleVar()
        self.volume = DoubleVar()
        self.f = DoubleVar()

        self.master = master

        self.master.geometry('1080x600+100+200')
        self.master.title("ещё одно окошко =)")
        self.label1 = Label(self.master, text="Введите параметры", fg='red').grid(row=0, column=2)

        self.nlabel = Label(self.master, text='количество каскадов').grid(row=1, column=0)
        self.deltaUnlabel = Label(self.master, text='приращение напряжение на N-ом каскаде').grid(row=2, column=0)
        self.Unlabel = Label(self.master, text='исходное напряжение').grid(row=3, column=0)
        self.Mlabel = Label(self.master, text='исходный параметр умножителя').grid(row=4, column=0)
        self.sigmalabel = Label(self.master, text='сигма').grid(row=5, column=0)
        self.wavelabel = Label(self.master, text='длина волны').grid(row=6, column=0)
        self.Tlabel = Label(self.master, text='температура').grid(row=7, column=0)
        self.receiver_templabel = Label(self.master, text='температура приемника').grid(row=8, column=0)
        self.absorption_coefflabel = Label(self.master, text='коэф-т поглощения').grid(row=9, column=0)
        self.frequency_bandlabel = Label(self.master, text='частота приема').grid(row=10, column=0)
        self.photodetector_arealabel = Label(self.master, text='площадь приемника').grid(row=11, column=0)
        self.ΔF_backlabel = Label(self.master, text='фоновая частота').grid(row=12, column=0)
        self.S_intReceiverRadiationlabel = Label(self.master, text='инт. чувств-ь').grid(row=13, column=0)
        self.ΔF_receiverRadiationlabel = Label(self.master, text='частота поглощения').grid(row=14, column=0)
        self.S_vlabel = Label(self.master, text='макс. инт. чувст-ь(рассчитывать)').grid(row=15, column=0)
        self.R0label = Label(self.master, text='вн. сопротивление').grid(row=16, column=0)
        self.Ilabel = Label(self.master, text='шум тока(рассчитывать)').grid(row=17, column=0)
        self.RLoadlabel = Label(self.master, text='сопр. нагрузки').grid(row=18, column=0)
        self.darkCurlabel = Label(self.master, text='темн. ток').grid(row=19, column=0)
        self.V_supplylabel = Label(self.master, text='напр. питания').grid(row=20, column=0)
        self.R_darklabel = Label(self.master, text='темн. сопр.').grid(row=21, column=0)
        self.tau_carrierlabel = Label(self.master, text='время жизни носителей').grid(row=22, column=0)
        self.nCarlabel = Label(self.master, text='кол-во носителей').grid(row=23, column=0)
        self.volumelabel = Label(self.master, text='шум фотослоя').grid(row=24, column=0)
        self.flabel = Label(self.master, text='частота модуляции').grid(row=25, column=0)


        self.nEntry = Entry(self.master, textvar=self.n).grid(row=1, column=2)
        self.deltaUnEntry = Entry(self.master, textvar=self.deltaUn).grid(row=2, column=2)
        self.UnEntry = Entry(self.master, textvar=self.Un).grid(row=3, column=2)
        self.MEntry = Entry(self.master, textvar=self.M).grid(row=4, column=2)
        self.sigmaEntry = Entry(self.master, textvar=self.sigma).grid(row=5, column=2)
        self.waveEntry = Entry(self.master, textvar=self.wave).grid(row=6, column=2)
        self.TEntry = Entry(self.master, textvar=self.T).grid(row=7, column=2)
        self.receiver_tempEntry = Entry(self.master, textvar=self.receiver_temp).grid(row=8, column=2)
        self.absorption_coeffEntry = Entry(self.master, textvar=self.absorption_coeff).grid(row=9, column=2)
        self.frequency_bandEntry = Entry(self.master, textvar=self.frequency_band).grid(row=10, column=2)
        self.photodetector_areaEntry = Entry(self.master, textvar=self.photodetector_area).grid(row=11, column=2)
        self.ΔF_backEntry = Entry(self.master, textvar=self.ΔF_back).grid(row=12, column=2)
        self.S_intReceiverRadiationEntry = Entry(self.master,
                                                 textvar=self.S_intReceiverRadiation).grid(row=13, column=2)
        self.ΔF_receiverRadiationEntry = Entry(self.master, textvar=self.ΔF_receiverRadiation).grid(row=14, column=2)
        self.S_vnEntry = Entry(self.master, textvar=self.S_v).grid(row=15, column=2)
        self.R0Entry = Entry(self.master, textvar=self.R0).grid(row=16, column=2)
        self.IEntry = Entry(self.master, textvar=self.I).grid(row=17, column=2)
        self.RLoadEntry = Entry(self.master, textvar=self.RLoad).grid(row=18, column=2)
        self.darkCurEntry = Entry(self.master, textvar=self.darkCur).grid(row=19, column=2)
        self.V_supplyEntry = Entry(self.master, textvar=self.V_supply).grid(row=20, column=2)
        self.R_darkEntry = Entry(self.master, textvar=self.R_dark).grid(row=21, column=2)
        self.tau_carrierEntry = Entry(self.master, textvar=self.tau_carrier).grid(row=22, column=2)
        self.nCarEntry = Entry(self.master, textvar=self.nCar).grid(row=23, column=2)
        self.volumeEntry = Entry(self.master, textvar=self.volume).grid(row=24, column=2)
        self.fEntry = Entry(self.master, textvar=self.f).grid(row=25, column=2)


        self.button1 = Button(self.master, text='зав. изм. коэф. ус-я', fg='blue',
                              command=self.deltaM).grid(row=17, column=5)
        self.button2 = Button(self.master, text='макс. инт. зав-ь', fg='blue',
                              command=self.Sv).grid(row=18, column=5)
        self.button3 = Button(self.master, text='флуктуация излучения', fg='blue',
                              command=self.fluct).grid(row=19, column=5)
        self.button4 = Button(self.master, text='дисперс.напр. рад. шума', fg='blue',
                              command=self.vRadSq).grid(row=20, column=5)
        self.button5 = Button(self.master, text='терм. шум напр.', fg='blue',
                              command=self.vol).grid(row=21, column=5)
        self.button6 = Button(self.master, text='терм. шум тока', fg='blue',
                              command=self.cur).grid(row=22, column=5)
        self.button7 = Button(self.master, text='дроб. шум тока', fg='blue',
                              command=self.shot).grid(row=23, column=5)
        self.button8 = Button(self.master, text='дроб. напр.', fg='blue', command=self.Vshot).grid(row=24, column=5)
        self.button9 = Button(self.master, text='шум лав. тока', fg='blue', command=self.avNoise).grid(row=24, column=5)
        self.button10 = Button(self.master, text='фон. шум напр.', fg='blue', command=self.Vgr).grid(row=25, column=5)
        self.button0 = Button(self.master, text='Назад', fg='blue', command=self.quit).grid(row=5, column=3)

    #V_gr = voltage_gr_noise(5, R_load, 100, 25, 1e-6, 1e18, 1e-9, 1000, frequency_band)
    #print('Фоновый шум напряжения: ', V_gr)

    def Vgr(self):
        global state
        M = self.M.get()
        darkCur = self.darkCur.get()
        frequency_band = self.frequency_band.get()

        I_shot = avalanche_shot_noise(darkCur, M, frequency_band)

        print(I_shot)
        state = "Шум лавинного тока =  " + str(I_shot)
        self.labelresult = Label(self.master, text=state).grid(row=9, column=10)
        return state

    def avNoise(self):
        global state
        M = self.M.get()
        darkCur = self.darkCur.get()
        frequency_band = self.frequency_band.get()

        I_shot = avalanche_shot_noise(darkCur, M, frequency_band)

        print(I_shot)
        state = "Шум лавинного тока =  " + str(I_shot)
        self.labelresult = Label(self.master, text=state).grid(row=9, column=10)
        return state

    def Vshot(self):
        global state
        RLoad = self.RLoad.get()
        I = self.I.get()
        frequency_band = self.frequency_band.get()

        V_shot = rms_current_shot_voltage(I, RLoad, frequency_band)

        print(V_shot)
        state = "Дробовое напряжение =  " + str(V_shot)
        self.labelresult = Label(self.master, text=state).grid(row=8, column=10)
        return state

    def shot(self):
        global state
        I = self.I.get()
        frequency_band = self.frequency_band.get()

        ShotNoise = rms_current_shot_noise(I, frequency_band)

        print(ShotNoise)
        state = "Дробовый шум тока =  " + str(ShotNoise)
        self.labelresult = Label(self.master, text=state).grid(row=7, column=10)
        return state

    def cur(self):
        global state
        I = self.I.get()
        frequency_band = self.frequency_band.get()

        ShotNoise = rms_current_shot_noise(I, frequency_band)

        print(ShotNoise)
        state = "Дробовый шум тока =  " + str(ShotNoise)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=10)
        return state

    def cur(self):
        global state
        T = self.T.get()
        R0 = self.R0.get()
        frequency_band = self.frequency_band.get()

        current = calc_thermal_noise_current(T, R0, frequency_band)

        print(current)
        state = "Термический шум тока =  " + str(current)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=10)
        return state

    def vol(self):
        global state
        T = self.T.get()
        R0 = self.R0.get()
        frequency_band = self.frequency_band.get()

        voltage2 = calc_thermal_noise_voltage(T, R0, frequency_band)

        print(voltage2)
        state = "Термический шум напряжения =  " + str(voltage2)
        self.labelresult = Label(self.master, text=state).grid(row=5, column=10)
        return state

    def vRadSq(self):
        global state
        ΔF_back = self.ΔF_back.get()
        S_intReceiverRadiation = self.S_intReceiverRadiation.get()
        ΔF_receiverRadiation = self.ΔF_receiverRadiation.get()
        S_v = self.S_v.get()

        V_radiation_squared = calculate_radiation_noise_variance(S_v, ΔF_back,
                                                                 S_intReceiverRadiation, ΔF_receiverRadiation)

        print(V_radiation_squared)
        state = "дисперсия напряжения рад. шума =  " + str(V_radiation_squared)
        self.labelresult = Label(self.master, text=state).grid(row=4, column=10)
        return state

    def fluct(self):
        global state
        receiver_temp = self.receiver_temp.get()
        absorption_coeff = self.absorption_coeff.get()
        frequency_band = self.frequency_band.get()
        photodetector_area = self.photodetector_area.get()
        T = self.T.get()

        result = noise_power(T, receiver_temp, absorption_coeff, frequency_band, photodetector_area)

        print(result)
        state = "флуктуация потоков излучения фона и премника =  " + str(result)
        self.labelresult = Label(self.master, text=state).grid(row=3, column=10)
        return state

    def Sv(self):
        global state

        def calculate_Un(wavelength):
            h = 6.62607015e-34  # Planck constant in J*s
            c = 299792458  # Speed of light in m/s
            q = 1.602176634e-19  # Elementary charge in C

            return (h * c / q) / wavelength

        sigma = self.sigma.get()
        Un = calculate_Un(self.wave.get())
        T = self.T.get()

        S_v = max_integral_sensitivity(sigma, Un, T)

        print(S_v)
        state = "Максимальная интегральная чувствительность =  " + str(S_v)
        self.labelresult = Label(self.master, text=state).grid(row=2, column=10)
        return state

    def deltaM(self):
        global state
        n = self.n.get()
        delta_U_n = self.deltaUn.get()
        U_n = self.Un.get()
        M = self.M.get()

        delta_M = delta_gain(n, delta_U_n, U_n, M)

        print(delta_M)
        state = "зависимость изменения коэф. усиления от изм. напр. питания =  " + str(delta_M)
        self.labelresult = Label(self.master, text=state).grid(row=1, column=10)
        return state

    def quit(self):
        self.master.destroy()