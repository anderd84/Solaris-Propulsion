from dataclasses import dataclass
from types import NoneType
import numpy as np
import matplotlib.pyplot as plt
import matrix_viewer as mv

from Cooling.material import DomainMaterial
from Cooling import material
from Nozzle import plots

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
        plt.ion()
        plt.show()

        for i in range(self.hpoints):
            for j in range(self.vpoints):

                if prevPercent < int(i * j / (self.hpoints * self.vpoints) * 200):
                    prevPercent = int(i * j / (self.hpoints * self.vpoints) * 200)
                    print(f"Progress: {prevPercent/2}%")
                    fig.axes[0].clear()
                    plots.PlotPlug(fig, plug)
                    plots.PlotPlug(fig, cowl)
                    self.ShowMaterialPlot(fig)
                    fig.canvas.draw()
                    fig.canvas.flush_events()

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
                



def SolveDomain(domain: Domain, tol: float = 1e-3):
    for i in range(domain.hpoints):
        for j in range(domain.vpoints):
            if i == 0 or j == 0 or i == domain.hpoints-1 or j == domain.vpoints-1:
                print("edge")
