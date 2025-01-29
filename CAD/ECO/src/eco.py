import tkinter as tk
from tkinter import ttk
import os
import yaml

def main():
    win = tk.Tk()
    win.title("ECO Form")

    win.geometry("800x600")

    winStyle = ttk.Style(win)
    win.tk.call('source', os.path.join(os.path.dirname(__file__),'theme/arc.tcl'))
    winStyle.theme_use('arc')
    win.configure(background="#FFFFFF")

    frame = ttk.Frame(win)
    frame.pack()

    label = ttk.Label(frame, text="ECO Form").pack()
    ttk.Label(frame, text="Name").pack()
    nameTextBox = ttk.Entry(frame).pack()


    win.mainloop()

def getPrevECOnumber():
    pass

















if __name__ == '__main__':
    main()