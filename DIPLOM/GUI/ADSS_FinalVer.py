from tkinter import *
import CalcGUI

# ЦЕНТРАЛЬНОЕ ОКНО
class MainPage:
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('СППР "Пиротех"')
        self.label1 = Label(self.master, text='Выберите интересующую вас функцию',
                            fg='red').grid(row=0, column=2)
        self.button1 = Button(self.master, text='расчет базовых параметров', fg='blue', bg="white",
                              bd=10, command=self.Selector).grid(row=6, column=2)
        self.button2 = Button(self.master, text='расчет спектральных параметров', fg='blue', bg="white", bd=10,
                              command=self.SpectreSelector).grid(row=8, column=2)
        self.button2 = Button(self.master, text='База данных фотоприемников', fg='blue', bg="white", bd=10,
                              command=self.DBPR).grid(row=10, column=2)
        self.button1 = Button(self.master, text='выход', fg='blue', bg="white", bd=10,
                              command=self.finish).grid(row=6, column=5)

    def Selector(self):
        root2 = Toplevel(self.master)
        myGui = CalcSelector2(root2)

    def SpectreSelector(self):
        root2 = Toplevel(self.master)
        myGui = Spectral(root2)

    def DBPR(self):
        root2 = Toplevel(self.master)
        myGui = DB(root2)

    def finish(self):
        self.master.destroy()


class Spectral(MainPage):
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('окошко)')


class DB(MainPage):
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('окошко)')

class CalcSelector2(MainPage):
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('окошко)')
        self.label1 = Label(self.master, text='Выберите, что хотите рассчитать',
                            fg='red').grid(row=0, column=2)
        self.button1 = Button(self.master, text='Фоторезистор', fg='blue',
                              command=self.PhotoResistor).grid(row=1, column=1)
        self.button2 = Button(self.master, text='Фотодиод', fg='blue', command=self.PhotoDiode).grid(row=2, column=1)
        self.button3 = Button(self.master, text='Фотоумножитель', fg='blue',
                              command=self.PhotoMultiplyer).grid(row=3, column=1)

    def PhotoResistor(self):
        root3 = Toplevel(self.master)
        myGui = CalcGUI.PhoRes(root3)

    def PhotoDiode(self):
        root3 = Toplevel(self.master)
        myGui = CalcGUI.PhoDio(root3)

    def PhotoMultiplyer(self):
        root3 = Toplevel(self.master)
        myGui = CalcGUI.PhoMul(root3)

    def finish(self):
        self.master.destroy()



class PlotSelector(MainPage):
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('для графиков)')
        self.label1 = Label(self.master, text='Выберите, что хотите рассчитать',
                            fg='red').grid(row=0, column=2)
        self.button1 = Button(self.master, text='Фототок', fg='blue', bg="white", bd=10,
                              command=self.PhotoCur).grid(row=1, column=1)
        self.button2 = Button(self.master, text='Ср. чувст-ь', fg='blue', bg="white", bd=10,
                              command=self.PhotoSens).grid(row=2, column=1)
        self.button3 = Button(self.master, text='Дифф. чувст-ь', fg='blue', bg="white", bd=10,
                              command=self.PhotoDifSens).grid(row=3, column=1)
        self.button3 = Button(self.master, text='Темн. ток', fg='blue', bg="white", bd=10,
                              command=self.DarkCur).grid(row=4, column=1)

    def DarkCur(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoDarkCur(root3)

    def PhotoDifSens(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoDiffSensivity(root3)

    def PhotoSens(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoSensivity(root3)

    def PhotoCur(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoCurrent(root3)

    def finish(self):
        self.master.destroy()


class CalcSelector(MainPage):
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('окошко)')
        self.label1 = Label(self.master, text='Выберите, что хотите рассчитать',
                            fg='red').grid(row=0, column=2)
        self.button1 = Button(self.master, text='Фототок', fg='blue', command=self.PhotoCur).grid(row=1, column=1)
        self.button2 = Button(self.master, text='Ср. чувст-ь', fg='blue',
                              command=self.PhotoSens).grid(row=2, column=1)
        self.button3 = Button(self.master, text='Дифф. чувст-ь', fg='blue',
                              command=self.PhotoDifSens).grid(row=3, column=1)
        self.button3 = Button(self.master, text='Темн. ток', fg='blue',
                              command=self.DarkCur).grid(row=4, column=1)

    def DarkCur(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoDarkCur(root3)

    def PhotoDifSens(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoDiffSensivity(root3)

    def PhotoSens(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoSensivity(root3)

    def PhotoCur(self):
        root3 = Toplevel(self.master)
        #myGui = PhotoCurrent(root3)

    def finish(self):
        self.master.destroy()



def main():
    root = Tk()
    myGuiWelcome = MainPage(root)
    root.mainloop()


if __name__ == '__main__':
    main()