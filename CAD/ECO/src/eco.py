import tkinter as tk
from tkinter import ttk
import os
from subprocess import run
import yaml
import git

config = yaml.load(open(os.path.join(os.path.dirname(__file__), 'config.yaml'), 'r'), Loader=yaml.FullLoader)

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

    getChangedFiles()

    # win.mainloop()

def getPrevECOnumber():
    pass



def getChangedFiles() -> set[str]:
    """
    Get the list of files that have been changed in the current commit that are tracked by ECOs

    Returns:
        set: A set of file paths that have been changed
    """
    repo = git.Repo(os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

    diffStaged = [diff.a_path for diff in repo.index.diff(repo.head.commit)]
    diffs = [diff.a_path for diff in repo.index.diff(None)]
    diffs = set(diffStaged) | set(diffs)

    cadChanges = set()

    for diff in diffs:
        _, ext = os.path.splitext(diff)
        print(ext)
        if ext.lower() in config["trackedFileExtensions"]:
            cadChanges.add(diff)

    return cadChanges



if __name__ == '__main__':
    main()