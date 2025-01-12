import numpy as np
import pint
from enum import Enum
from typing import Any
import tempfile
import random
import string
from alive_progress import alive_bar, alive_it

from cooling.domain import DomainMC, DomainPoint, MaterialType
from general.units import unitReg, Q_

class DomainMMAP(DomainMC):
    attributes: list = []
    memmaps: dict

    units: dict

    workingFolder = tempfile.mkdtemp()
    prefix = str(''.join(random.choices(string.ascii_letters, k=4)))

    def __init__(self, domain: DomainMC):
        self.x0 = domain.x0
        self.r0 = domain.r0
        self.width = domain.width
        self.height = domain.height
        self.hpoints = domain.hpoints
        self.vpoints = domain.vpoints
        self.xstep = domain.xstep
        self.rstep = domain.rstep

        print("Loading domain")
        self.attributes = list(DomainPoint(0, 0, 0).__dict__.keys())
        self.memmaps = {}
        self.units = {}

        for attr in alive_it(self.attributes):
            print(f"Transferring {attr}")
            testAttr = domain.array[0,0].__getattribute__(attr)
            dtype = 'float64'
            if isinstance(testAttr, pint.Quantity):
                self.units[attr] = str(testAttr.units)
                t = [[domain.array[i,j].__getattribute__(attr).to(testAttr.units).magnitude for j in range(domain.hpoints)] for i in range(domain.vpoints)]
            elif isinstance(testAttr, Enum):
                t = [[domain.array[i,j].__getattribute__(attr).value for j in range(domain.hpoints)] for i in range(domain.vpoints)]
                dtype = 'int'
            elif isinstance(testAttr, tuple):
                t = [[[a for a in domain.array[i,j].__getattribute__(attr)] for j in range(domain.hpoints)] for i in range(domain.vpoints)]
                dtype = 'int'
            else:
                t = [[domain.array[i,j].__getattribute__(attr) for j in range(domain.hpoints)] for i in range(domain.vpoints)]
            shape = (domain.vpoints, domain.hpoints) if not isinstance(testAttr, tuple) else (domain.vpoints, domain.hpoints, len(testAttr))
            self.memmaps[attr] = np.memmap(f'{self.workingFolder}/{self.prefix}_{attr}.dat', dtype=dtype, mode='w+', shape=shape)
            self.memmaps[attr][:] = t[:]
            del t

    def __getattribute__(self, name: str) -> Any:
        try:
            return super().__getattribute__(name)
        except AttributeError:
            if name in self.attributes:
                if name in self.units:
                    return Q_(self.memmaps[name][:], self.units[name])
                return self.memmaps[name]
            else:
                return super().__getattribute__(name)
            
    def __setattr__(self, name: str, value: Any) -> None:
        if name in ['attributes', 'memmaps', 'units', 'workingFolder', 'prefix'] or name in super().__dir__():
            super().__setattr__(name, value)
        else:
            if name in self.attributes:
                if isinstance(value, pint.Quantity):
                    self.memmaps[name][:] = value.to(Q_(self.units[name])).magnitude
                elif isinstance(value, Enum):
                    self.memmaps[name][:] = value.value
                else:
                    self.memmaps[name][:] = value
                self.memmaps[name].flush()
            else:
                super().__setattr__(name, value)

    def toDomain(self):
        domain = DomainMC(self.x0, self.r0, self.width, self.height, self.xstep)
        with alive_bar(self.vpoints*self.hpoints, title="Setting point data") as bar:
            for i in range(self.vpoints):
                for j in range(self.hpoints):
                    for attr in self.attributes:
                        if attr in self.units:
                            domain.array[i,j].__setattr__(attr, Q_(self.memmaps[attr][i,j], self.units[attr]))
                        elif isinstance(self.memmaps[attr][i,j], np.memmap):
                            domain.array[i,j].__setattr__(attr, tuple(self.memmaps[attr][i,j]))
                        else:
                            domain.array[i,j].__setattr__(attr, self.memmaps[attr][i,j])
                    bar()
        return domain
    
    def setMEM(self, row, col, name, value):
        if isinstance(value, pint.Quantity):
            self.memmaps[name][row, col] = value.to(Q_(self.units[name])).magnitude
        elif isinstance(value, Enum):
            self.memmaps[name][row, col] = value.value
        else:
            self.memmaps[name][row, col] = value
        self.memmaps[name].flush()

class SparseDomain(DomainMC):
    attributes: list = []
    points: dict

    units: dict

    def __init__(self, domain: DomainMC, poi: tuple[int, int]):
        self.x0 = domain.x0
        self.r0 = domain.r0
        self.width = domain.width
        self.height = domain.height
        self.hpoints = domain.hpoints
        self.vpoints = domain.vpoints
        self.xstep = domain.xstep
        self.rstep = domain.rstep

        self.attributes = list(DomainPoint(0, 0, 0).__dict__.keys())
        self.points = {}
        self.units = {}
        self.poi = poi

        prevPoi = domain.array[poi].previousFlow

        pointList = [poi, (poi[0] - 1, poi[1]), (poi[0] + 1, poi[1]), (poi[0], poi[1] - 1), (poi[0], poi[1] + 1),
                     prevPoi, (prevPoi[0] - 1, prevPoi[1]), (prevPoi[0] + 1, prevPoi[1]), (prevPoi[0], prevPoi[1] - 1), (prevPoi[0], prevPoi[1] + 1)]

        for attr in self.attributes:
            testAttr = domain.array[0,0].__getattribute__(attr)
            # if isinstance(testAttr, pint.Quantity):
            #     self.units[attr] = str(testAttr.units)
            self.points[attr] = {}
            for p in pointList:
                p = (max(min(p[0], domain.vpoints - 1), 0), max(min(p[1], domain.hpoints - 1), 0))
                if isinstance(testAttr, pint.Quantity):
                    t = domain.array[p].__getattribute__(attr)#.to(testAttr.units).magnitude
                elif isinstance(testAttr, Enum):
                    t = domain.array[p].__getattribute__(attr).value
                elif isinstance(testAttr, tuple):
                    t = [a for a in domain.array[p].__getattribute__(attr)]
                else:
                    t = domain.array[p].__getattribute__(attr)
                self.points[attr][p] = t

    def __getattribute__(self, name: str) -> Any:
        try:
            return super().__getattribute__(name)
        except AttributeError:
            if name in self.attributes:
                if name in self.units:
                    return Q_(self.points[name], self.units[name])
                return self.points[name]
            else:
                return super().__getattribute__(name)
            
    def __setattr__(self, name: str, value: Any) -> None:
        if name in ['attributes', 'points', 'units', 'poi'] or name in super().__dir__():
            super().__setattr__(name, value)
        else:
            if name in self.attributes:
                if isinstance(value, pint.Quantity):
                    self.points[name] = value.to(Q_(self.units[name])).magnitude
                elif isinstance(value, Enum):
                    self.points[name] = value.value
                else:
                    self.points[name][:] = value
            else:
                super().__setattr__(name, value)
    
    def refreshUnits(self):
        for key in self.units.keys():
            if isinstance(self.points[key], pint.Quantity):
                self.points[key] = Q_(self.points[key].magnitude, self.units[key])