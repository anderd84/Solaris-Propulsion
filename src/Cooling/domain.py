from dataclasses import dataclass
from tracemalloc import start
from types import NoneType
import concurrent
from fluids.gas import Gas
import numpy as np
np.product = np.prod
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
from General.units import Q_, unitReg

mcQ = Queue()
SHAREDMEMNAME = 'CoolingDomain'

@dataclass
class DomainPoint:
    x: float
    r: float
    area: float
    material: DomainMaterial = DomainMaterial.FREE
    temperature: Q_ = Q_(70, unitReg.degF)
    velocity: Q_ = Q_(0, unitReg.ft/unitReg.s)

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
        
        self.array = np.empty((vpoints, hpoints), dtype=DomainPoint)

        self.x0 = x0
        self.r0 = r0
        self.width = width
        self.height = height
        self.hpoints = hpoints
        self.vpoints = vpoints
        self.xstep = width/(hpoints-1)
        self.rstep = height/(vpoints-1)
        for i in range(vpoints):
            for j in range(hpoints):
                self.array[i,j] = DomainPoint(x0 + j*self.xstep, r0 - i*self.rstep, self.xstep*self.rstep)

    def DefineMaterials(self, cowl: np.ndarray, coolant: np.ndarray, chamber: np.ndarray, plug: np.ndarray, max_cores = mp.cpu_count() - 1):
        MAX_CORES = max_cores
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
            for i in range(self.vpoints):
                for j in range(self.hpoints):
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
        matarr = np.array([[point.temperature.magnitude for point in row] for row in self.array])

        # extent = [xarr[0,0]-self.xstep, xarr[-1,-1]+self.xstep, rarr[-1,-1]-self.rstep, rarr[0,0]+self.rstep]
        extent = [xarr[0,0]-self.xstep/2, xarr[-1,-1]+self.xstep/2, rarr[-1,-1]-self.rstep/2, rarr[0,0]+self.rstep/2]
        ax = fig.axes[0]
        # contf = ax.contourf(xarr, rarr, matarr, levels=[0, 1, 4] , colors=['white', 'blue', 'red'])
        # contf = ax.contourf(xarr, rarr, matarr)
        ax.imshow(matarr, extent=extent, origin='upper', cmap='jet')
        xcells = np.linspace(self.x0 - self.xstep/2, self.x0 + self.width + self.xstep/2, self.hpoints+1)
        rcells = np.linspace(self.r0 + self.rstep/2, self.r0 - self.height - self.rstep/2, self.vpoints+1)
        xl, rl = np.meshgrid(xcells, rcells)
        ax.plot(xl, rl, 'k', linewidth=0.25)
        ax.plot(np.transpose(xl), np.transpose(rl), 'k', linewidth=0.25)

    def ShowStatePlot(self, fig: plt.Figure):
        xarr = np.array([[point.x for point in row] for row in self.array])
        rarr = np.array([[point.r for point in row] for row in self.array])
        matarr = np.array([[point.temperature.magnitude for point in row] for row in self.array])

        ax = fig.axes[0]
        contf = ax.contourf(xarr, rarr, matarr, cmap='jet')
        xcells = np.linspace(self.x0 - self.xstep/2, self.x0 + self.width + self.xstep/2, self.hpoints+1)
        rcells = np.linspace(self.r0 + self.rstep/2, self.r0 - self.height - self.rstep/2, self.vpoints+1)
        xl, rl = np.meshgrid(xcells, rcells)
        ax.plot(xl, rl, 'k', linewidth=0.25)
        ax.plot(np.transpose(xl), np.transpose(rl), 'k', linewidth=0.25)

    def AssignChamberTemps(self, chamber: np.ndarray, exhaust: Gas, startPoint: tuple, endPoint: tuple, chamberWallRadius: Q_):
        flowAngle = np.arctan((endPoint[1] - startPoint[1])/(endPoint[0] - startPoint[0]))
        phi = np.pi/2 - flowAngle

        farAway1 = (startPoint[0] - chamberWallRadius.magnitude*np.sin(flowAngle), startPoint[1] + chamberWallRadius.magnitude*np.cos(flowAngle))
        farAway2 = (startPoint[0] + chamberWallRadius.magnitude*np.sin(flowAngle), startPoint[1] - chamberWallRadius.magnitude*np.cos(flowAngle))
        startPointU = material.intersectPolyAt(chamber, startPoint, farAway1)
        startPointL = material.intersectPolyAt(chamber, startPoint, farAway2)



        originX = np.linspace(startPoint[0], endPoint[0], self.hpoints)
        originR = np.linspace(startPoint[1], endPoint[1], self.hpoints)

        iStart, jStart = self.ChamberStartCell()

        for j in range(jStart, self.hpoints):
            for i in range(iStart, self.vpoints):
                if self.array[i,j].material != DomainMaterial.CHAMBER:
                    self.array[i,j].material = DomainMaterial.COOLANT
                    break
                if self.lineInCell(startPointU, startPointL, i, j):
                    jStart = j if i == iStart else jStart
                    continue
                self.array[i,j].temperature = exhaust.stagTemp
                self.array[i,j].velocity = Q_(0, unitReg.foot/unitReg.sec)
            if self.lineInCell(startPointU, startPointL, i, j):
                break

        self.array[iStart,jStart].material = DomainMaterial.COOLANT       

        for oX, oR in zip(originX, originR):
            farAway1 = (oX - chamberWallRadius.magnitude*np.sin(flowAngle), oR + chamberWallRadius.magnitude*np.cos(flowAngle))
            farAway2 = (oX + chamberWallRadius.magnitude*np.sin(flowAngle), oR - chamberWallRadius.magnitude*np.cos(flowAngle))
            startPointU = material.intersectPolyAt(chamber, (oX, oR), farAway1)
            startPointL = material.intersectPolyAt(chamber, (oX, oR), farAway2)
            # plt.plot([startPointU[0], startPointL[0]], [startPointU[1], startPointL[1]], 'rx-')
            area = np.pi/np.sin(phi) * (startPointU[1]**2 - startPointL[1]**2)
            
    def isInCell(self, point, row, col):
        xmin = self.array[row, col].x - self.xstep/2
        xmax = self.array[row, col].x + self.xstep/2
        ymin = self.array[row, col].r - self.rstep/2
        ymax = self.array[row, col].r + self.rstep/2
        return xmin <= point[0] and point[0] <= xmax and ymin <= point[1] and point[1] <= ymax
        
    def lineInCell(self, point1, point2, row, col):
        linex = np.linspace(point1[0], point2[0], self.hpoints)
        liner = np.linspace(point1[1], point2[1], self.hpoints)
        line = np.array([[linex[i], liner[i]] for i in range(self.hpoints)])
        for point in line:
            if self.isInCell(point, row, col):
                return True

    def ChamberStartCell(self):
        for i in range(self.vpoints):
            for j in range(self.hpoints):
                if self.array[i,j].material == DomainMaterial.CHAMBER:
                    return (i, j)
        

    def DumpFile(self, filename):
        with open(filename, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

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
