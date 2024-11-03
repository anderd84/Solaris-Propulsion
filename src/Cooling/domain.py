from dataclasses import dataclass
from types import NoneType
import concurrent
import numpy as np
import matplotlib.pyplot as plt
import matrix_viewer as mv

import multiprocessing as mp
from multiprocessing import shared_memory
from multiprocessing import Queue
import concurrent.futures
from time import sleep
import pickle

from icecream import ic

from Cooling.material import DomainMaterial
from Cooling import material
from Nozzle import plots

mcQ = Queue()
SHAREDMEMNAME = 'CoolingDomain'

@dataclass
class DomainPoint:
    x: float
    r: float
    material: DomainMaterial
    temperature: float
    area: float

@dataclass
class Domain:
    array: np.ndarray
    x0: float
    r0: float
    width: float
    height: float
    xstep: float
    rstep: float
    hpoints: int
    vpoints: int

    def __init__(self, x0, r0, width, height, dx=.1, dy=.1):
        hpoints = int(width/dx) + 1
        vpoints = int(height/dy) + 1
        self.array = np.empty((hpoints, vpoints), dtype=DomainPoint)
        self.x0 = x0
        self.r0 = r0
        self.width = width
        self.height = height
        self.hpoints = hpoints
        self.vpoints = vpoints
        self.xstep = width/(hpoints-1)
        self.rstep = height/(vpoints-1)
        for i in range(hpoints):
            for j in range(vpoints):
                self.array[i][j] = DomainPoint(x0 + i*self.xstep, r0 - j*self.rstep, DomainMaterial.FREE, 0, self.xstep*self.rstep)

    def DefineMaterials(self, cowl: np.ndarray, coolant: np.ndarray, chamber: np.ndarray, plug: np.ndarray, fig: plt.Figure):
        prevPercent = 0
        # plt.ion()
        # plt.show()

        for i in range(self.hpoints):
            for j in range(self.vpoints):

                if prevPercent < int(i * j / (self.hpoints * self.vpoints) * 200):
                    prevPercent = int(i * j / (self.hpoints * self.vpoints) * 200)
                    print(f"Progress: {prevPercent/2}%")
                    # fig.axes[0].clear()
                    # plots.PlotPlug(fig, plug)
                    # plots.PlotPlug(fig, cowl)
                    # self.ShowMaterialPlot(fig)
                    # fig.canvas.draw()
                    # fig.canvas.flush_events()

                if material.isIntersect(self.array[i][j], coolant, (self.width, self.height)):
                    self.array[i,j].material = DomainMaterial.COOLANT
                    continue
                if material.isIntersect(self.array[i][j], cowl, (self.width, self.height)):
                    self.array[i,j].material = DomainMaterial.COWL
                    continue
                if material.isIntersect(self.array[i][j], chamber, (self.width, self.height)):
                    self.array[i,j].material = DomainMaterial.CHAMBER
                    continue
                if material.isIntersect(self.array[i][j], plug, (self.width, self.height)):
                    self.array[i,j].material = DomainMaterial.PLUG
                    continue


    def ShowMaterialPlot(self, fig: plt.Figure):
        xarr = np.array([[point.x for point in row] for row in self.array])
        rarr = np.array([[point.r for point in row] for row in self.array])
        matarr = np.transpose(np.array([[point.material.value for point in row] for row in self.array]))

        ax = fig.axes[0]
        # contf = ax.contourf(xarr, rarr, matarr, levels=[0, 1, 4] , colors=['white', 'blue', 'red'])
        # contf = ax.contourf(xarr, rarr, matarr)
        ax.imshow(matarr, extent=(self.x0, self.x0 + self.width, self.r0 - self.height, self.r0), origin='upper', cmap='jet')
        xcells = np.linspace(self.x0 - self.xstep/2, self.x0 + self.width + self.xstep/2, self.hpoints+1)
        rcells = np.linspace(self.r0 + self.rstep/2, self.r0 - self.height - self.rstep/2, self.vpoints+1)
        xl, rl = np.meshgrid(xcells, rcells)
        # ax.plot(xl, rl, 'k', linewidth=0.25)
        # ax.plot(np.transpose(xl), np.transpose(rl), 'k', linewidth=0.25)
@dataclass
class DomainMC:
    array: np.ndarray
    x0: float
    r0: float
    width: float
    height: float
    xstep: float
    rstep: float
    hpoints: int
    vpoints: int

    def __init__(self, x0, r0, width, height, ds = .1):        
        hpoints = int(width/ds) + 1
        vpoints = int(height/ds) + 1
        
        self.array = np.empty((hpoints, vpoints), dtype=DomainPoint)

        self.x0 = x0
        self.r0 = r0
        self.width = width
        self.height = height
        self.hpoints = hpoints
        self.vpoints = vpoints
        self.xstep = width/(hpoints-1)
        self.rstep = height/(vpoints-1)
        for i in range(hpoints):
            for j in range(vpoints):
                self.array[i][j] = DomainPoint(x0 + i*self.xstep, r0 - j*self.rstep, DomainMaterial.FREE, 0, self.xstep*self.rstep)

    def DefineMaterials(self, cowl: np.ndarray, coolant: np.ndarray, chamber: np.ndarray, plug: np.ndarray, fig: plt.Figure):
        MAX_CORES = 1
        shm = None
        try:
            shm = shared_memory.SharedMemory(create=True, size=self.array.nbytes, name=SHAREDMEMNAME)
        except FileExistsError:
            shm = shared_memory.SharedMemory(name=SHAREDMEMNAME)
            shm.unlink()
            shm = shared_memory.SharedMemory(create=True, size=self.array.nbytes, name=SHAREDMEMNAME)
        
        newarray = np.ndarray(self.array.shape, dtype=DomainPoint, buffer=shm.buf)
        newarray[:] = self.array[:]

        self.array = newarray.copy()

        # use a pool of processes to parallelize the computation
        with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_CORES) as executor:
            q = []
            for i in range(self.hpoints):
                for j in range(self.vpoints):
                    # q.append((i, j))
                    mcQ.put((i, j))

            futures = []
            print(f"Starting parallel computation with {MAX_CORES} cores")
            jump = int(len(q)/MAX_CORES)
            for i in range(MAX_CORES):
                futures.append(executor.submit(EvalProcess, i, shm, self.array.shape, (self.width, self.height), coolant, cowl, chamber, plug))
                # futures.append(executor.submit(EvalProcess, MAX_CORES - i - 1, q[i*jump:(i+1)*jump], self.shm, self.array.shape, (self.width, self.height), coolant, cowl, chamber, plug))
            print("Tasks submitted")
            for future in concurrent.futures.as_completed(futures):
                res = future.result()
                for i, j, mat in res:
                    self.array[i,j].material = mat

        print("Parallel computation done")
        shm.close()
        shm.unlink()

    def ShowMaterialPlot(self, fig: plt.Figure):
        xarr = np.array([[point.x for point in row] for row in self.array])
        rarr = np.array([[point.r for point in row] for row in self.array])
        matarr = np.transpose(np.array([[point.material.value for point in row] for row in self.array]))

        ax = fig.axes[0]
        # contf = ax.contourf(xarr, rarr, matarr, levels=[0, 1, 4] , colors=['white', 'blue', 'red'])
        # contf = ax.contourf(xarr, rarr, matarr)
        ax.imshow(matarr, extent=(self.x0, self.x0 + self.width, self.r0 - self.height, self.r0), origin='upper', cmap='jet')
        xcells = np.linspace(self.x0 - self.xstep/2, self.x0 + self.width + self.xstep/2, self.hpoints+1)
        rcells = np.linspace(self.r0 + self.rstep/2, self.r0 - self.height - self.rstep/2, self.vpoints+1)
        xl, rl = np.meshgrid(xcells, rcells)
        # ax.plot(xl, rl, 'k', linewidth=0.25)
        # ax.plot(np.transpose(xl), np.transpose(rl), 'k', linewidth=0.25)

    def DumpFile(self, filename):
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def LoadFile(filename):
        with open(filename, "rb") as f:
            return pickle.load(f)

def EvalProcess(pn, shm, shape, size, coolant, cowl, chamber, plug):
    res = []
    domain = np.ndarray(shape, dtype=DomainPoint, buffer=shm.buf)
    print(f"Starting process {pn + 1}")
    prevPercent = 0
    total = np.prod(shape)
    while mcQ.qsize() > 0:
        i, j = mcQ.get()
        res.append(AssignMaterial(domain, i, j, size, coolant, cowl, chamber, plug))
        if pn == 0:
            curPercent = int((total - mcQ.qsize())/total * 100)
            if prevPercent < curPercent:
                prevPercent = curPercent
                print(f"Progress: {prevPercent}%")
    return res

def AssignMaterial(domain, i, j, size, coolant, cowl, chamber, plug):
    if material.isIntersect(domain[i][j], coolant, size):
        return (i, j, DomainMaterial.COOLANT)
    if material.isIntersect(domain[i][j], cowl, size):
        return (i, j, DomainMaterial.COWL)
    if material.isIntersect(domain[i][j], chamber, size):
        return (i, j, DomainMaterial.CHAMBER)
    if material.isIntersect(domain[i][j], plug, size):
        return (i, j, DomainMaterial.PLUG)
    return (i, j, DomainMaterial.FREE)
    