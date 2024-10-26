from dataclasses import dataclass
from types import NoneType
import numpy as np

from Cooling.material import DomainMaterial
from Cooling import material

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

    def __init__(self, x0, r0, width, height, hpoints: int = 100, vpoints: int = 100):
        self.array = np.empty((hpoints, vpoints), dtype=DomainPoint)
        self.x0 = x0
        self.r0 = r0
        self.width = width
        self.height = height
        self.hpoints = hpoints
        self.vpoints = vpoints
        self.xstep = width/hpoints
        self.rstep = height/vpoints
        for i in range(hpoints):
            for j in range(vpoints):
                self.array[i][j] = DomainPoint(x0 + i*self.xstep, r0 - j*self.rstep, DomainMaterial.FREE, 0, self.xstep*self.rstep)

    def DefineMaterials(self, cowl: np.ndarray, coolant: np.ndarray, chamber: np.ndarray, plug: np.ndarray):
        for i in range(self.hpoints):
            for j in range(self.vpoints):
                if material.isIntersect(self.array[i][j], cowl, (self.width, self.height)):
                    self.array[i][j].material = DomainMaterial.COWL
                if material.isIntersect(self.array[i][j], coolant, (self.width, self.height)):
                    self.array[i][j].material = DomainMaterial.COOLANT
                if material.isIntersect(self.array[i][j], chamber, (self.width, self.height)):
                    self.array[i][j].material = DomainMaterial.CHAMBER
                if material.isIntersect(self.array[i][j], plug, (self.width, self.height)):
                    self.array[i][j].material = DomainMaterial.PLUG
                
                



def SolveDomain(domain: Domain, tol: float = 1e-3):
    for i in range(domain.hpoints):
        for j in range(domain.vpoints):
            if i == 0 or j == 0 or i == domain.hpoints-1 or j == domain.vpoints-1:
                print("edge")
