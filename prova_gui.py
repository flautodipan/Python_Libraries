#%%

import tkinter as tk

def first_print():
    text = 'Porca la Madonna!'
    text_output = tk.Label(window, text = text)
    text_output.grid(row=0, column = 1)


window = tk.Tk()
window.geometry('600x500')
window.title('PorcoDio!')
window.configure(background = 'red')
first_button = tk.Button(text = 'Bestemmia!', command = first_print)
first_button.grid(row = 0, column = 0, )

if __name__ == '__main__':
    window.mainloop()

# %%
