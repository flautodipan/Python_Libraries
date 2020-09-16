#%%

import tkinter as tk


def first_print():
    text = 'Porca la Madonna!'
    text_output = tk.Label(window, text = text)
    #text_output.pack(side = tk.TOP, fill = tk.BOTH, expand = False)
    text_output.grid(row=0, column = 1)
def second_print():
    text_output = tk.Label(window, text = "Hai scelto di bestemmiare con: " + in_put.get() )
    text_output.grid(row=1, column = 1)



window = tk.Tk()
in_put = tk.Entry(window, width = 20)
window.iconbitmap('/home/tesista/Desktop/preview.jpg')
in_put.grid(row = 1, column = 0)
window.geometry('600x500')
window.title('PorcoDio!')
#window.configure()
first_button = tk.Button(text = 'Bestemmia!', command = first_print, background = 'red')
second_button = tk.Button(text = 'Lancia la tua bestemmia', command = second_print, background= 'green')
first_button.grid(row = 0, column = 0, )
second_button.grid(row = 2, column = 0)

if __name__ == '__main__':
    window.mainloop()

# %%
