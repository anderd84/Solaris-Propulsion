import datetime
import tkinter as tk
from tkinter import ttk
import os
from datetime import datetime
import yaml
import git
from glob import glob
from datetime import datetime

config = yaml.load(open(os.path.join(os.path.dirname(__file__), 'config.yaml'), 'r'), Loader=yaml.FullLoader)
topLevelDir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
repo = git.Repo(topLevelDir)

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

    ttk.Label(frame, text="ECO Form").grid(row=0, column=0, columnspan=2)
    ttk.Label(frame, text="Name").grid(row=1, column=0)
    nameTextBox = ttk.Entry(frame)
    nameTextBox.grid(row=1, column=1)
    ttk.Label(frame, text="Date").grid(row=2, column=0)
    dateTextBox = ttk.Entry(frame)
    dateTextBox.grid(row=2, column=1)
    dateTextBox.insert(0, datetime.now().strftime("%m/%d/%Y"))

    

    

    win.mainloop()

def getPrevECOnumber():
    cadDirectory = os.path.join(topLevelDir, "CAD")
    result = [y for x in os.walk(cadDirectory) for y in glob(os.path.join(x[0], 'ECO_*.md'))]
    return max([int(file[-6:-3]) for file in result])

def getChangedFiles() -> set[str]:
    """
    Get the list of files that have been changed in the current commit that are tracked by ECOs

    Returns:
        set: A set of file paths that have been changed
    """
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

def getModifiedDate(file: str) -> str:
    return datetime.fromtimestamp(os.stat(file).st_mtime).date()

if __name__ == '__main__':
    main()