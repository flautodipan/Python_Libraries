#%%
# 
from tkinter import *
from PIL import ImageTk, Image

root = Tk()
root.title('Images')

my_img = ImageTk.PhotoImage(Image.open('../elenina.png'))
my_label = Label(image = my_img)
my_label.pack()


exit_button = Button(text = 'Exit program', command = root.quit)
exit_button.pack()

if __name__ == '__main__':
    root.mainloop()

# %%
