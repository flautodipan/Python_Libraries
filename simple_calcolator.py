#%%


from tkinter import *
import numpy as np

window = Tk()
window.title('Simple calculator')

in_put = Entry(window, width = 35, borderwidth = 5)
in_put.grid(row = 0, column = 0, columnspan = 3, padx = 10, pady = 10)
in_put.insert(0, 'Calcoliamo!')

def button_click(number):
    if in_put.get() == 'Calcoliamo!': in_put.delete(0, END)
    current = in_put.get()
    in_put.delete(0, END)
    in_put.insert(0, str(current) + str(number))

def button_clear():
    in_put.delete(0, END)

def button_add():
    global prev
    global operation

    prev = int(in_put.get())
    operation = '+'

    print('sommo')
    in_put.delete(0, END)


def button_multiply():
    global prev
    global operation
    prev = int(in_put.get())
    operation = '*'
    print('moltiplico')
    in_put.delete(0, END)

def button_equal():
    succ = int(in_put.get())
    in_put.delete(0,END)

    if operation =='+':
        in_put.insert(0, str(prev+succ))
    elif operation == '*':
        in_put.insert(0, str(prev*succ))
    






buttons = []
operations_buttons = []

buttons.append(Button(window, text = '0', padx= 40, pady = 20, command = lambda: button_click(0)))
buttons[0].grid(row = 1, column = 1)
buttons.append(Button(window, text = '1', padx= 40, pady = 20, command = lambda: button_click(1)))
buttons.append(Button(window, text = '2', padx= 40, pady = 20, command = lambda: button_click(2)))
buttons.append(Button(window, text = '3', padx= 40, pady = 20, command = lambda: button_click(3)))
buttons.append(Button(window, text = '4', padx= 40, pady = 20, command = lambda: button_click(4)))
buttons.append(Button(window, text = '5', padx= 40, pady = 20, command = lambda: button_click(5)))
buttons.append(Button(window, text = '6', padx= 40, pady = 20, command = lambda: button_click(6)))
buttons.append(Button(window, text = '7', padx= 40, pady = 20, command = lambda: button_click(7)))
buttons.append(Button(window, text = '8', padx= 40, pady = 20, command = lambda: button_click(8)))
buttons.append(Button(window, text = '9', padx= 40, pady = 20, command = lambda: button_click(9)))

operations = ['C', '+', '*', '=', ]
fun_operations = [button_clear, button_add, button_multiply, button_equal]



for ii, (kk,jj) in zip(range(0, 10, 1), np.ndindex(3,3)):
    buttons[ii+1].grid(row = kk + 2 , column = jj )

for ll, op, c in zip(range(len(operations)), operations, fun_operations):

    operations_buttons.append(Button(window, text = op, padx= 39, pady = 20, command = c))
    operations_buttons[ll].grid(row = ll+1, column = 3)

if __name__ == '__main__':
    window.mainloop()

# %%
