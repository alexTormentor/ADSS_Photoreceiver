import tkinter as tk
from tkinter.constants import DISABLED, NORMAL
from collections import namedtuple
from tests import *

User = namedtuple("User", ["username", "password", "user_type"])



class PhotoCur(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.geometry("720x420")
        self.C = tk.Label(self, text="постоянная, специфичная для фоторезистора",
                     font="impact 15")
        self.U = tk.Label(self, text="напряжение на фоторезисторе (в вольтах)",
                     font="impact 15")
        self.W = tk.Label(self, text="энергия", font="impact 15")
        self.Time = tk.Label(self, text="время", font="impact 15")
        self.ConvCur = tk.Label(self, text="Ток, преобразованный из светового потока",
                           font="impact 15")
        self.C.place(x=10, y=10)
        self.U.place(x=10, y=50)
        self.W.place(x=10, y=90)
        self.Time.place(x=10, y=130)
        self.ConvCur.place(x=10, y=170)
        self.CValue = tk.StringVar()
        self.UValue = tk.StringVar()
        self.WValue = tk.StringVar()
        self.TimeValue = tk.StringVar()
        self.CEntry = tk.Entry(self, textvariable=self.CValue, font="impact 15", width=15)
        self.UEntry = tk.Entry(self, textvariable=self.UValue, font="impact 15", width=15)
        self.WEntry = tk.Entry(self, textvariable=self.WValue, font="impact 15", width=15)
        self.TimeEntry = tk.Entry(self, textvariable=self.TimeValue, font="impact 15", width=15)
        self.CEntry.place(x=450, y=10)
        self.UEntry.place(x=450, y=50)
        self.WEntry.place(x=450, y=90)
        self.TimeEntry.place(x=450, y=130)

    def aboba(self):
        self.ABOBA = tk.Button(text="Рассчитать", font="impact 15", bg="white",
               command=self.Calculation, bd=10).place(x=45, y=30)

    def Calculation(self):
        C = self.CEntry.get()
        C1 = float(C)
        U = self.UEntry.get()
        W = self.WEntry.get()
        Time = self.TimeEntry.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        current = photoresistor_current(C1, U, a, F, b)
        tk.Label(text=f"{current}", font="impact 15").place(x=450, y=170)


class CalcForm(tk.Toplevel):
    def __init__(self, parent, calc_type):
        super().__init__(parent)
        self.calc_type = calc_type
        self.geometry("720x420")

        self.geometry("720x420")
        self.C = tk.Label(self, text="постоянная, специфичная для фоторезистора",
                          font="impact 15")
        self.U = tk.Label(self, text="напряжение на фоторезисторе (в вольтах)",
                          font="impact 15")
        self.W = tk.Label(self, text="энергия", font="impact 15")
        self.Time = tk.Label(self, text="время", font="impact 15")
        self.ConvCur = tk.Label(self, text="Ток, преобразованный из светового потока",
                                font="impact 15")
        self.C.place(x=10, y=10)
        self.U.place(x=10, y=50)
        self.W.place(x=10, y=90)
        self.Time.place(x=10, y=130)
        self.ConvCur.place(x=10, y=170)
        self.CValue = tk.StringVar()
        self.UValue = tk.StringVar()
        self.WValue = tk.StringVar()
        self.TimeValue = tk.StringVar()
        self.CEntry = tk.Entry(self, textvariable=self.CValue, font="impact 15", width=15)
        self.UEntry = tk.Entry(self, textvariable=self.UValue, font="impact 15", width=15)
        self.WEntry = tk.Entry(self, textvariable=self.WValue, font="impact 15", width=15)
        self.TimeEntry = tk.Entry(self, textvariable=self.TimeValue, font="impact 15", width=15)
        self.CEntry.place(x=450, y=10)
        self.UEntry.place(x=450, y=50)
        self.WEntry.place(x=450, y=90)
        self.TimeEntry.place(x=450, y=130)
        btn = tk.Button(self, text="Фототок", bg="white",
                        command=self.open2, bd=10)
        btn.grid(row=3, columnspan=2)

    def open2(self):
        C = self.CEntry.get()
        C1 = float(C)
        U = float(self.UEntry.get())
        W = float(self.WEntry.get())
        Time = float(self.TimeEntry.get())
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        self.current2 = photoresistor_current(C1, U, a, F, b)
        tk.Label(text=f"{self.current2}", font="impact 15").place(x=450, y=170)


class App2(tk.Tk):
    def __init__(self):
        super().__init__()
        calc_types = ("Расчёт физических параметров", "В РАЗРАБОТКЕ")
        self.calc_type = tk.StringVar()
        self.calc_type.set(calc_types[0])


        label = tk.Label(self, text="Выберите интересующую функцию")
        radios = [tk.Radiobutton(self, text=t, value=t,
                                 variable=self.calc_type) for t in calc_types]
        btn = tk.Button(self, text="Выбрать", command=self.open_window)

        label.pack(padx=10, pady=10)
        for radio in radios:
            radio.pack(padx=10, anchor=tk.W)
        btn.pack(pady=10)

    def open_window(self):
        window = CalcForm(self, self.calc_type.get())
        calcs = window.open2()
        print(calcs)



if __name__ == "__main__":
    app = App2()
    app.mainloop()