from tkinter import *
from tests import *


class Welcome:
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('Welcome')
        self.label1 = Label(self.master, text='welcome to the wages calculator GUI', fg='red').grid(row=0, column=2)
        self.button1 = Button(self.master, text='ok', fg='blue', command=self.gotowages).grid(row=6, column=2)
        self.button1 = Button(self.master, text='quit', fg='blue', command=self.finish).grid(row=6, column=3)

    def gotowages(self):
        root2 = Toplevel(self.master)
        myGui = Wages(root2)

    def finish(self):
        self.master.destroy()


state = ''


class Wages(Welcome):
    def __init__(self, master):
        # Welcome.__init__(self,*args,**kwargs)
        self.mhours = DoubleVar()
        self.salaryh = DoubleVar()
        self.master = master
        # self.salary=salary
        self.master.geometry('500x200+100+200')
        self.master.title("Wages Calculator")
        self.label1 = Label(self.master, text="welcome to salary calculator", fg='red').grid(row=0, column=2)
        self.label2 = Label(self.master, text='Enetr your salary per hour').grid(row=3, column=0)
        self.label3 = Label(self.master, text='Enter number of hours worked').grid(row=4, column=0)
        self.mysalary = Entry(self.master, textvar=self.salaryh).grid(row=3, column=2)
        self.myhours = Entry(self.master, textvar=self.mhours).grid(row=4, column=2)
        self.button1 = Button(self.master, text='calculatesalary', fg='blue', command=self.calculatesalary).grid(row=5,
                                                                                                                 column=2)
        self.button2 = Button(self.master, text='Back', fg='blue', command=self.quit).grid(row=5, column=3)
        self.button3 = Button(self.master, text='page3', fg='blue', command=self.page3).grid(row=5, column=4)

    def calculatesalary(self):
        global state
        hours = self.mhours.get()
        # print (hours)
        hsal = self.salaryh.get()
        salary = hours * hsal
        print(salary)
        state = "your salary is " + str(salary)
        self.labelresult = Label(self.master, text=state).grid(row=7, column=2)
        # return state

    def page3(self):
        root3 = Toplevel(self.master)
        myPage = Page(root3)
        print(state)

    def quit(self):
        self.master.destroy()


# ---------------------------
class Page(Wages):
    def __init__(self, master):
        self.equation = StringVar()
        self.master = master
        self.master.geometry('500x200+100+200')
        self.master.title("trial page")
        self.label1 = Label(self.master, text="bibhu", fg='red').grid(row=0, column=2)
        self.label2 = Label(self.master, text="bibhu prasanna behera", fg='red').grid(row=0, column=3)
        self.button3 = Button(self.master, text='back', fg='blue', command=self.quit).grid(row=5, column=4)
        self.entry1 = Entry(self.master, textvar=self.equation, fg='red').grid(row=3, column=3)
        self.equation.set(state)
        self.button4 = Button(self.master, text='display', fg='blue', command=self.display).grid(row=7, column=4)

    def quit(self):
        self.master.destroy()

    def display(self):
        self.label3 = Label(self.master, text=self.equation.get(), fg='blue').grid(row=9, column=3)



class MainPage:
    def __init__(self, master):
        self.master = master
        self.master.geometry('400x200+100+200')
        self.master.title('СППР "Пиротех"')
        self.label1 = Label(self.master, text='Выберите интересующую вас функцию', fg='red').grid(row=0, column=2)
        self.button1 = Button(self.master, text='расчет базовых параметров', fg='blue', command=self.Basic).grid(row=6, column=2)
        self.button1 = Button(self.master, text='выход', fg='blue', command=self.finish).grid(row=6, column=3)

    def Basic(self):
        root2 = Toplevel(self.master)
        myGui = BasicCalc(root2)

    def finish(self):
        self.master.destroy()


class BasicCalc(MainPage):
    def __init__(self, master):
        # Welcome.__init__(self,*args,**kwargs)
        self.C = DoubleVar()
        self.U = DoubleVar()
        self.W = DoubleVar()
        self.Time = DoubleVar()

        self.mhours = DoubleVar()
        self.salaryh = DoubleVar()
        self.master = master
        # self.salary=salary

        self.master.geometry('500x200+100+200')
        self.master.title("AAAAAAAA")
        self.label1 = Label(self.master, text="Введите параметры", fg='red').grid(row=0, column=2)
        self.Clabel = Label(self.master, text='Введите константу').grid(row=1, column=0)
        self.Ulabel = Label(self.master, text='Введите напряжение').grid(row=2, column=0)
        self.Wlabel = Label(self.master, text='Введите энергию').grid(row=3, column=0)
        self.Timelabel = Label(self.master, text='Введите время').grid(row=4, column=0)
        self.CEntry = Entry(self.master, textvar=self.C).grid(row=1, column=2)
        self.UEntry = Entry(self.master, textvar=self.U).grid(row=2, column=2)
        self.WEntry = Entry(self.master, textvar=self.W).grid(row=3, column=2)
        self.TimeEntry = Entry(self.master, textvar=self.Time).grid(row=4, column=2)



        self.button1 = Button(self.master, text='рассчитать фототок', fg='blue', command=self.calculatesalary).grid(row=9,
                                                                                                                 column=2)
        self.button2 = Button(self.master, text='Back', fg='blue', command=self.quit).grid(row=5, column=3)
        self.button3 = Button(self.master, text='page3', fg='blue', command=self.page3).grid(row=5, column=4)

    def calculatesalary(self):
        global state
        C = self.C.get()
        U = self.U.get()
        W = self.W.get()
        Time = self.Time.get()
        voltage = [0, 0.5, 1, 1.5, 2, 2.5]
        current = [0, 0.12, 0.3, 0.4, 0.5, 0.6]
        intensity = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
        a = 0
        a = calculate_nonlinearity_coefficient(voltage, current)
        F = calculate_radiation_flux(energy=W, time=Time)
        b = calculate_light_nonlinearity_coefficient(intensity, current)
        current2 = photoresistor_current(C, U, a, F, b)
        print(current2)
        state = "Ваш фототок =  " + str(current2)
        self.labelresult = Label(self.master, text=state).grid(row=6, column=2)
        # return state

    def page3(self):
        root3 = Toplevel(self.master)
        myPage = Page(root3)
        print(state)

    def quit(self):
        self.master.destroy()



def main():
    root = Tk()
    myGuiWelcome = MainPage(root)
    root.mainloop()


if __name__ == '__main__':
    main()

