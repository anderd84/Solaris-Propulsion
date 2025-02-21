import datetime
import tkinter as tk
from tkinter import ttk
import os
from datetime import datetime
import yaml
import git
from glob import glob
from datetime import datetime
from enum import StrEnum

class ECOEffect(StrEnum):
    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"
    F = "F"

config = yaml.load(open(os.path.join(os.path.dirname(__file__), 'config.yaml'), 'r'), Loader=yaml.FullLoader)
topLevelDir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
repo = git.Repo(topLevelDir)

class ChangeTable(ttk.Frame):
    effect: list[ttk.Combobox] = []
    use: list[ttk.Checkbutton] = []
    data: dict[str, list] = {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data["names"] = list(getChangedFiles())
        self.data["oldRev"] = [getPrevECOnumber() for _ in self.data["names"]]
        self.data["newRev"] = [getPrevECOnumber() + 1 for _ in self.data["names"]]
        self.data["effect"] = [tk.StringVar(value=ECOEffect.A) for _ in self.data["names"]]
        self.data["date"] = [getModifiedDate(file) for file in self.data["names"]]
        self.data["use"] = [tk.BooleanVar(value=True) for _ in self.data["names"]]

        self._BuildTable()

    def _BuildTable(self):
        tableHeight = len(self.data["names"])
        ttk.Label(self, text="Include").grid(row=0, column=0)
        ttk.Label(self, text="Name").grid(row=0, column=2)
        ttk.Label(self, text="Old Rev").grid(row=0, column=4)
        ttk.Label(self, text="New Rev").grid(row=0, column=6)
        ttk.Label(self, text="Effect").grid(row=0, column=8)
        ttk.Label(self, text="Date").grid(row=0, column=10)
        for i in range(1,10,2):
            ttk.Separator(self, orient="vertical").grid(row=0, column=i, rowspan=tableHeight+2, sticky='ns')

        ttk.Separator(self, orient="horizontal").grid(row=1, column=0, columnspan=11, sticky='ew')

        self.use = [None for _ in range(tableHeight)]
        self.effect = [None for _ in range(tableHeight)]
        self.date = [None for _ in range(tableHeight)]

        for i in range(tableHeight):
            self.use[i] = ttk.Checkbutton(self, variable=self.data["use"][i])
            self.use[i].grid(row=i+2, column=0)
            ttk.Label(self, text=self.data["names"][i]).grid(row=i+2, column=2)
            ttk.Label(self, text=self.data["oldRev"][i]).grid(row=i+2, column=4)
            ttk.Label(self, text=self.data["newRev"][i]).grid(row=i+2, column=6)
            self.effect[i] = ttk.Combobox(self, values=[effect.value for effect in ECOEffect], textvariable=self.data["effect"][i])
            self.effect[i].grid(row=i+2, column=8)
            self.effect[i].current(0)
            ttk.Label(self, text=self.data["date"][i]).grid(row=i+2, column=10)





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

    ttk.Label(frame, text="ECO Form").grid(row=0, column=0, columnspan=10)
    ttk.Label(frame, text="Name").grid(row=1, column=0)
    nameTextBox = ttk.Entry(frame)
    nameTextBox.grid(row=1, column=1)
    ttk.Label(frame, text="Date").grid(row=2, column=0)
    dateTextBox = ttk.Entry(frame)
    dateTextBox.grid(row=2, column=1)
    dateTextBox.insert(0, datetime.now().strftime("%m/%d/%Y"))

    ttk.Separator(frame, orient="horizontal").grid(row=3, column=0, columnspan=10, sticky='ew')

    ttk.Label(frame, text="Information on Effects").grid(row=4, column=0, columnspan=10)
    ttk.Label(frame, text=
              "Please select the appropriate effectivity level of changes:\n- A: Clerical – No impact to production\n- B: No Delay – Incorporate as routine parts availability allows\n- C: Immediate – All new production must be built to new revision\n- D: Rework Stock – All production inventory must be reworked to the new revision\n- E: Rework Finished Goods – No shipment to old revision\n- F: Scrap Stock – Scrap all materials produced to the old revision"
              ).grid(row=5, column=0, columnspan=10)

    ttk.Separator(frame, orient="horizontal").grid(row=6, column=0, columnspan=10, sticky='ew')

    changeTable = ChangeTable(frame)
    changeTable.grid(row=7, column=0, columnspan=10)
    

    

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