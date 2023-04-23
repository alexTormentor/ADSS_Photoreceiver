import tkinter as tk
from tests import *

class Demo1:
    def __init__(self, master):
        self.master = master
        self.master.geometry("400x400")
        self.frame = tk.Frame(self.master)
        self.butnew("Window 1", "ONE", Demo2)
        self.butnew("Window 2", "TWO", Demo3)
        self.frame.pack()

    def butnew(self, text, number, _class):
        tk.Button(self.frame, text = text, width = 25, command = lambda: self.new_window(number, _class)).pack()

    def new_window(self, number, _class):
        self.newWindow = tk.Toplevel(self.master)
        _class(self.newWindow, number)


class Demo2:
    def __init__(self, master, number):
        self.master = master
        self.master.geometry("400x400+400+400")
        self.frame = tk.Frame(self.master)
        self.quitButton = tk.Button(self.frame, text = 'Quit', width = 25, command = self.close_windows)
        self.label = tk.Label(master, text=f"this is window number {number}")
        self.label.pack()
        self.quitButton.pack()
        self.frame.pack()

    def close_windows(self):
        self.master.destroy()

class Demo3:
    def __init__(self, master, number):
        self.master = master
        self.master.geometry("400x400+400+400")
        self.frame = tk.Frame(self.master)
        self.quitButton = tk.Button(self.frame, text = 'Quit', width = 25, command = self.close_windows)
        self.label = tk.Label(master, text=f"this is window number {number}")
        self.label.pack()
        self.label2 = tk.Label(master, text="THIS IS HERE TO DIFFERENTIATE THIS WINDOW")
        self.label2.pack()
        self.quitButton.pack()
        self.frame.pack()

    def close_windows(self):
        self.master.destroy()



class BasicCalc:
    def __init__(self, master, number):
        self.master = master
        self.master.geometry("400x400+400+400")
        self.frame = tk.Frame(self.master)
        self.quitButton = tk.Button(self.frame, text='Закрыть', width=25, command=self.close_windows)
        self.butnew("Фототок", "базовые параметры", Calculator)
        self.label = tk.Label(master, text=f"Здесь вы можете рассчитать {number} фотоприбора")
        self.label.pack()
        self.quitButton.pack()
        self.frame.pack()

    def butnew(self, text, number, _class):
        tk.Button(self.frame, text = text, bg="white", bd=10, width = 25, command = lambda: self.new_window(number, _class)).pack()

    def close_windows(self):
        self.master.destroy()

    def new_window(self, number, _class):
        self.newWindow = tk.Toplevel(self.master)
        _class(self.newWindow, number)

class Calculator:
    def __init__(self, master, number):
        self.master = master
        self.master.geometry("720x420")
        self.frame = tk.Frame(self.master)
        self.C = tk.Label(master, text="постоянная, специфичная для фоторезистора",
                          font="impact 15")
        self.U = tk.Label(master, text="напряжение на фоторезисторе (в вольтах)",
                          font="impact 15")
        self.W = tk.Label(master, text="энергия", font="impact 15")
        self.Time = tk.Label(master, text="время", font="impact 15")
        self.ConvCur = tk.Label(master, text="Ток, преобразованный из светового потока",
                                font="impact 15")
        self.C.pack(padx=0, pady=1)
        self.U.pack(padx=0, pady=3)
        self.W.pack(padx=0, pady=4)
        self.Time.pack(padx=0, pady=5)
        self.ConvCur.pack(padx=0, pady=6)

        self.CValue = tk.StringVar()
        self.UValue = tk.StringVar()
        self.WValue = tk.StringVar()
        self.TimeValue = tk.StringVar()
        self.CEntry = tk.Entry(master, textvariable=self.CValue, font="impact 15", width=15)
        self.UEntry = tk.Entry(master, textvariable=self.UValue, font="impact 15", width=15)
        self.WEntry = tk.Entry(master, textvariable=self.WValue, font="impact 15", width=15)
        self.TimeEntry = tk.Entry(master, textvariable=self.TimeValue, font="impact 15", width=15)
        self.CEntry.pack(padx=4, pady=1)
        self.UEntry.pack(padx=4, pady=3)
        self.WEntry.pack(padx=4, pady=4)
        self.TimeEntry.pack(padx=4, pady=5)

        self.quitButton = tk.Button(self.frame, text='Закрыть', width=25, command=self.close_windows)
        self.butnew("Фототок", C=self.CEntry, U=self.UEntry, W=self.WEntry, Time=self.TimeEntry, _class=Calculator)

        self.frame.pack()

    def butnew(self, text, C, U, W, Time, _class):
        tk.Button(self.frame, text = text, C = C, U = U, W = W, Time=Time, bg="white", bd=10, width = 25,
                  command = lambda: self.calcs(_class=Calculator)).pack()

    def calcs(self, text, _class):
        C = float(self.CEntry.get())
        U = float(self.UEntry.get())
        W = float(self.WEntry.get())
        Time = float(self.TimeEntry.get())
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        current2 = photoresistor_current(C, U, a, F, b)
        tk.Label(text=f"{current2}", font="impact 15").pack(padx=0, pady=6)

    def close_windows(self):
        self.master.destroy()

class MainWindow:
    def __init__(self, master):
        self.master = master
        self.master.geometry("400x400")
        self.master.title('СППР "Пиротех"')
        self.frame = tk.Frame(self.master)
        self.butnew("Расчеты", "базовые параметры", BasicCalc)
        self.frame.pack()

    def butnew(self, text, number, _class):
        tk.Button(self.frame, text = text, bg="white", bd=10, width = 25, command = lambda: self.new_window(number, _class)).pack()


    def new_window(self, number, _class):
        self.newWindow = tk.Toplevel(self.master)
        _class(self.newWindow, number)

def main(): 
    root = tk.Tk()
    app = MainWindow(root)
    root.mainloop()

if __name__ == '__main__':
    main()