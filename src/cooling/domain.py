from dataclasses import dataclass
import shutil
import time
from nozzle.analysis import CharacteristicPoint
import numpy as np
import pint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import griddata
from typing import Any
import multiprocessing as mp
from icecream import ic
import joblib
from alive_progress import alive_bar, alive_it

from cooling.material import DomainMaterial, MaterialType, CoolantType
from cooling import material
from fluids import fluid, gas
from fluids.gas import Gas
from general.units import Q_, unitReg
from general import design

@dataclass
class CoolingChannel:
    upperContour: np.ndarray
    lowerContour: np.ndarray    

@dataclass
class CoolingLoop:
    landWidth: pint.Quantity
    numChannels: int
    mdot: pint.Quantity
    fluid: CoolantType

@dataclass
class DomainPoint:
    x: float
    r: float
    area: pint.Quantity = Q_(0, unitReg.inch**2)
    material: DomainMaterial = DomainMaterial.FREE
    border: bool = False
    temperature: pint.Quantity = Q_(70, unitReg.degF).to(unitReg.degR)
    pressure: pint.Quantity = Q_(12.1, unitReg.psi)
    velocity: pint.Quantity = Q_(0, unitReg.ft/unitReg.s)
    hydraulicDiameter: pint.Quantity = Q_(0, unitReg.inch)
    previousFlow: tuple[int, int] = (-1,-1)
    futureFlow: tuple[int, int] = (-1,-1)
    flowHeight: pint.Quantity = Q_(0, unitReg.inch)
    id: int = -1

    def getState(self, state: str):
        try:
            s = self.__getattribute__(state)
        except AttributeError:
            raise AttributeError(f"DomainPoint has no attribute {state}")
        if isinstance(s, pint.Quantity):
            return s.magnitude
        return s

class DomainMC:
    array: np.ndarray[DomainPoint]
    x0: float
    r0: float
    width: float
    height: float
    xstep: float
    rstep: float
    hpoints: int
    vpoints: int
    coolingLoops: dict[int, CoolingLoop]

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
        self.coolingLoops = {}
        print("Creating domain")
        with alive_bar(vpoints*hpoints) as bar:
            for i in range(vpoints):
                for j in range(hpoints):
                    self.array[i,j] = DomainPoint(x0 + j*self.xstep, r0 - i*self.rstep, Q_(self.xstep*self.rstep, unitReg.inch**2))
                    bar()
        print("Domain created")

    def DefineMaterials(self, cowl: np.ndarray, chamber: np.ndarray, plug: np.ndarray, max_cores = mp.cpu_count() - 1):
        MAX_CORES = max_cores
        print(f"using {MAX_CORES} cores")
        tic = time.perf_counter()

        print("Starting processes")
        with joblib.Parallel(n_jobs=MAX_CORES, verbose=100, return_as='generator') as parallel:

            with alive_bar(self.vpoints*self.hpoints) as bar:
                print("Starting processes")

                outputs = parallel(joblib.delayed(EvalMaterialProcess2)(i, self.hpoints, i, (self.x0, self.xstep, self.r0, self.rstep), (self.width, self.height), cowl, chamber, plug) for i in range(self.vpoints))
                print(outputs)

                print("Assinging outputs")
                for out in outputs:
                    for i, j, mat in out:
                        self.array[i,j].material = mat
                        bar()
            
            try:
                shutil.rmtree('./.work')
            except OSError as e:
                print("Error: %s - %s." % (e.filename, e.strerror))

        print("Parallel computation done")
        print("assigning borders")
        self.AssignBorders()
        toc = time.perf_counter()
        print("material defined")
        print(f"Time to define materials: {toc - tic}")

    def ContourPlot(self, fig: plt.Figure, attr: str, omitMaterials: list = []):
        xarr = np.array([[point.x for point in row] for row in self.array])
        rarr = np.array([[point.r for point in row] for row in self.array])
        
        matarr = np.array([[point.getState(attr) if point.material not in omitMaterials else np.nan for point in row] for row in self.array])

        ax = fig.axes[0]
        contf = ax.contourf(xarr, rarr, matarr, 100, cmap='jet')
        units = str(self.array[0,0].__getattribute__(attr).units) if isinstance(self.array[0,0].__getattribute__(attr), pint.Quantity) else "unitless"
        fig.colorbar(contf, ax=ax).set_label(f"{attr} ({units})")
        print("done!")

    def NodePlot(self, fig: plt.Figure, attr: str, omitMaterials: list = []):
        xarr = np.array([[point.x for point in row] for row in self.array])
        rarr = np.array([[point.r for point in row] for row in self.array])
        
        matarr = np.array([[point.getState(attr) if point.material not in omitMaterials else np.nan for point in row] for row in self.array])

        extent = [xarr[0,0]-self.xstep/2, xarr[-1,-1]+self.xstep/2, rarr[-1,-1]-self.rstep/2, rarr[0,0]+self.rstep/2]
        ax = fig.axes[0]
        im = ax.imshow(matarr, extent=extent, origin='upper', cmap='jet')
        units = str(self.array[0,0].__getattribute__(attr).units) if isinstance(self.array[0,0].__getattribute__(attr), pint.Quantity) else "unitless"
        fig.colorbar(im, ax=ax).set_label(f"{attr} ({units})")

    def RelationPlot(self, fig: plt.Figure):
        ax = fig.axes[0]
        for i in range(self.vpoints):
            for j in range(self.hpoints):
                flow = self.array[i,j].futureFlow
                if (flow[0] != -1 and flow[1] != -1) or (self.array[i,j].material == DomainMaterial.COOLANT_WALL):
                    # ax.quiver(self.array[i,j].x, self.array[i,j].r, self.array[flow].x - self.array[i,j].x, self.array[flow].r - self.array[i,j].r, scale=1, scale_units='xy', angles='xy', color='r', width=0.002)
                    ax.quiver(self.array[flow].x, self.array[flow].r, self.array[i,j].x - self.array[flow].x, self.array[i,j].r - self.array[flow].r, scale=1, scale_units='xy', angles='xy', color='r', width=0.002)
                flow = self.array[i,j].previousFlow
                if flow[0] != -1:
                    if self.array[i,j].material == DomainMaterial.COOLANT_BULK:
                        ax.quiver(self.array[flow].x, self.array[flow].r, self.array[i,j].x - self.array[flow].x, self.array[i,j].r - self.array[flow].r, scale=1, scale_units='xy', angles='xy', color='k', width=0.002)
                        # ax.quiver(self.array[i,j].x, self.array[i,j].r, self.array[i,j].x - self.array[flow].x, self.array[i,j].r - self.array[flow].r, scale=1, scale_units='xy', angles='xy', color='k', width=0.002)
                    else:
                        ax.quiver(self.array[flow].x, self.array[flow].r, self.array[i,j].x - self.array[flow].x, self.array[i,j].r - self.array[flow].r, scale=1, scale_units='xy', angles='xy', color='b', width=0.002)
                        # ax.quiver(self.array[i,j].x, self.array[i,j].r, self.array[i,j].x - self.array[flow].x, self.array[i,j].r - self.array[flow].r, scale=1, scale_units='xy', angles='xy', color='b', width=0.002)   


    def ShowCellPlot(self, fig: plt.Figure):
        ax = fig.axes[0]
        xcells = np.linspace(self.x0 - self.xstep/2, self.x0 + self.width + self.xstep/2, self.hpoints+1)
        rcells = np.linspace(self.r0 + self.rstep/2, self.r0 - self.height - self.rstep/2, self.vpoints+1)
        xl, rl = np.meshgrid(xcells, rcells)
        ax.plot(xl, rl, 'k', linewidth=0.25)
        ax.plot(np.transpose(xl), np.transpose(rl), 'k', linewidth=0.25)

    def NewCoolantLoop(self, landWidth: pint.Quantity, numChannels: int, mdot: pint.Quantity, fluid: CoolantType):
        loopID = len(self.coolingLoops)
        self.coolingLoops[loopID] = CoolingLoop(landWidth, numChannels, mdot, fluid)
        return loopID

    def AssignCoolantFlow(self, coolant: CoolingChannel, upperWall: bool, initialPressure: pint.Quantity, loopID: int):
        # first assign wall
        if loopID not in self.coolingLoops:
            raise ValueError(f"Loop ID {loopID} not defined")
        landWidth = self.coolingLoops[loopID].landWidth
        numChannels = self.coolingLoops[loopID].numChannels

        inputPoints = len(coolant.upperContour)
        previousWall = (-1,-1)
        previousFlow = (-1,-1)
        # areas type
        # initial wall assignment
        wallContour = coolant.upperContour if upperWall else coolant.lowerContour
        previousMat = self.array[self.CoordsToCell(wallContour[0].x, wallContour[0].r)].material
        inlet = (-1,-1)
        for i in alive_it(range(inputPoints - 1)):
            dist = np.sqrt((wallContour[i].x - wallContour[i+1].x)**2 + (wallContour[i].r - wallContour[i+1].r)**2)
            steps = max(int(dist/min(self.xstep, self.rstep) * 1.5), 5)
            xc = np.linspace(wallContour[i].x, wallContour[i+1].x, steps)[:-1]
            rc = np.linspace(wallContour[i].r, wallContour[i+1].r, steps)[:-1]

            for j in range(steps - 1):
                if xc[j] > self.x0 + self.width or xc[j] < self.x0:
                    break
                if rc[j] < self.r0 - self.height or rc[j] > self.r0:
                    break
                wallPoint = self.CoordsToCell(xc[j], rc[j])
                if (i==0 and j==0):
                    self.array[wallPoint].material = DomainMaterial.COOLANT_INLET
                    inlet = wallPoint
                else:
                    self.array[wallPoint].material = DomainMaterial.COOLANT_WALL if self.array[wallPoint].material != DomainMaterial.COOLANT_INLET else DomainMaterial.COOLANT_INLET
                
                if wallPoint[0] != previousWall[0] or wallPoint[1] != previousWall[1]:
                    previousFlow = previousWall
                self.array[wallPoint].previousFlow = previousFlow
                self.array[previousFlow].futureFlow = wallPoint if self.array[wallPoint].material != DomainMaterial.COOLANT_INLET else (-1,-1)
                previousWall = wallPoint
        # delete bad points
        point = inlet
        while self.array[point].futureFlow != point:
            # print(point)
            nextPoint = self.array[point].futureFlow
            offset = (-1, 0) if upperWall else (1, 0)
            if point[0] + offset[0] < 0 or point[0] + offset[0] >= self.vpoints or point[1] + offset[1] < 0 or point[1] + offset[1] >= self.hpoints:
                point = nextPoint
                continue 
            if self.array[point[0] + offset[0], point[1] + offset[1]].material == DomainMaterial.COOLANT_WALL and (self.array[point[0], point[1] - 1].material == DomainMaterial.COOLANT_WALL or self.array[point[0], point[1] + 1].material == DomainMaterial.COOLANT_WALL):
                self.array[self.array[point].futureFlow].previousFlow = self.array[point].previousFlow
                self.array[self.array[point].previousFlow].futureFlow = self.array[point].futureFlow
                self.array[point].previousFlow = (-1,-1)
                self.array[point].futureFlow = (-1,-1)
                self.array[point].material = previousMat
            point = nextPoint

        previousWall = (-1,-1)
        previousFlow = (-1,-1)
        wallPoint = (-1,-1)

        pointMap = {}

        for i in alive_it(range(inputPoints - 1)):
            dist1 = np.sqrt((coolant.lowerContour[i].x - coolant.lowerContour[i+1].x)**2 + (coolant.lowerContour[i].r - coolant.lowerContour[i+1].r)**2)
            dist2 = np.sqrt((coolant.upperContour[i].x - coolant.upperContour[i+1].x)**2 + (coolant.upperContour[i].r - coolant.upperContour[i+1].r)**2)
            dist = max(dist1, dist2)

            steps = max(int(dist/min(self.xstep, self.rstep) * 1.5), 5)

            xl = np.linspace(coolant.lowerContour[i].x, coolant.lowerContour[i+1].x, steps)[:-1]
            rl = np.linspace(coolant.lowerContour[i].r, coolant.lowerContour[i+1].r, steps)[:-1]

            xu = np.linspace(coolant.upperContour[i].x, coolant.upperContour[i+1].x, steps)[:-1]
            ru = np.linspace(coolant.upperContour[i].r, coolant.upperContour[i+1].r, steps)[:-1]

            # plt.plot([xl, xu], [rl, ru], '-r', linewidth=.25)

            for j in range(steps - 1):
                if xl[j] > self.x0 + self.width or xu[j] > self.x0 + self.width:
                    break
                if xl[j] < self.x0 or xu[j] < self.x0:
                    break

                if rl[j] < self.r0 - self.height or ru[j] < self.r0 - self.height:
                    break
                if rl[j] > self.r0 or ru[j] > self.r0:
                    break

                possibleWall = self.CoordsToCell(xu[j], ru[j]) if upperWall else self.CoordsToCell(xl[j], rl[j])
                if self.array[possibleWall].material == DomainMaterial.COOLANT_WALL or self.array[possibleWall].material == DomainMaterial.COOLANT_INLET:
                    wallPoint = possibleWall

                cells = self.cellsOnLine((xl[j], rl[j]), (xu[j], ru[j]))

                # plt.plot([xl[j], xu[j]], [rl[j], ru[j]], '-b', linewidth=.25)

                landSectorAngle = landWidth/Q_(rl[j], unitReg.inch)
                channelSectorAngle = (2*np.pi/numChannels) - landSectorAngle
                h = Q_(np.sqrt((xl[j] - xu[j])**2 + (rl[j] - ru[j])**2), unitReg.inch)
                totalArea = np.pi*h*Q_(ru[j] + rl[j], unitReg.inch)
                channelArea = totalArea*channelSectorAngle/(2*np.pi)
                perim = (Q_(rl[j], unitReg.inch) + Q_(ru[j], unitReg.inch))*channelSectorAngle + 2*h
                hydroD = 4*channelArea/perim

                for cellPos in cells:
                    if (i==0 and j==0):
                        self.array[cellPos].material = DomainMaterial.COOLANT_INLET

                    edgeCase = (xl[j] >= self.x0 + self.width - self.xstep or xu[j] >= self.x0 + self.width - self.xstep)
                    edgeCase = edgeCase or (xl[j] <= self.x0 + self.xstep or xu[j] <= self.x0 + self.xstep)
                    edgeCase = edgeCase or (rl[j] <= self.r0 - self.height + self.rstep or ru[j] <= self.r0 - self.height + self.rstep)
                    edgeCase = edgeCase or (rl[j] >= self.r0 - self.rstep or ru[j] >= self.r0 - self.rstep)
                    if (i==inputPoints-2 and j==steps-2) or edgeCase:
                        self.array[cellPos].material = DomainMaterial.COOLANT_OUTLET

                    if self.array[cellPos].material not in MaterialType.COOLANT | {DomainMaterial.COOLANT_INLET, DomainMaterial.COOLANT_OUTLET}:
                        self.array[cellPos].material = DomainMaterial.COOLANT_BULK
                        self.array[cellPos].previousFlow = wallPoint
                        pointMap.setdefault(wallPoint, set()).add(cellPos)


                    self.array[cellPos].pressure = initialPressure
                    self.array[cellPos].flowHeight = h
                    self.array[cellPos].hydraulicDiameter = hydroD
                    self.array[cellPos].area = channelArea
                    self.array[cellPos].id = loopID
                
        print("done")
        print("assigning borders")
        self.AssignBorders(loopID)
        # adjust non wall side
        for wallPoint in pointMap:
            borderCells = [cell for cell in pointMap[wallPoint] if self.array[cell].border and self.array[cell].material == DomainMaterial.COOLANT_BULK]
            for i, cell in enumerate(borderCells[1:]):
                self.array[cell].futureFlow = borderCells[i] # i is the index behind current position cause enumerate starts at 0
                self.array[cell].border = False
        
        print("done")

    def GuessChannelState(self, loopID: int, endTemp: pint.Quantity):
        # start at the left side, scroll down and right until notice a channel wall material
        # follow the previous flow chain until channel inlet material
        startPoint = (-1,-1)
        for j in range(self.hpoints):
            for i in range(self.vpoints):
                if self.array[i,j].id == loopID and self.array[i,j].material == DomainMaterial.COOLANT_WALL:
                    startPoint = (i,j)
                    break
            if startPoint != (-1,-1):
                break
        
        if self.array[startPoint].material != DomainMaterial.COOLANT_WALL:
            raise ValueError(f"No coolant found for id {loopID}")
        
        while self.array[startPoint].material != DomainMaterial.COOLANT_INLET:
            startPoint = self.array[startPoint].previousFlow
            if startPoint == (-1,-1):
                raise ValueError(f"No coolant inlet found for id {loopID}")
        
        # now we have the inlet point, future flow can be used to traverse the entire coolant loop

        currentPoint = startPoint
        deltaT = endTemp - self.array[startPoint].temperature
        pointData = []
        totalScore = 0
        channelPoints = []
        while self.array[currentPoint].futureFlow != currentPoint: # go until terminates
            channelPoints.append(currentPoint)
            currentPoint = self.array[currentPoint].futureFlow

        for currentPoint in alive_it(channelPoints):
            # score will be based on the max adjacent temperature 
            offsets = [(-1, 0), (0, 1), (1, 0), (0, -1)]
            maxTempSolid = Q_(0, unitReg.degR)

            for o in offsets:
                if currentPoint[0] + o[0] < 0 or currentPoint[0] + o[0] >= self.vpoints:
                    continue
                if currentPoint[1] + o[1] < 0 or currentPoint[1] + o[1] >= self.hpoints:
                    continue
                t = self.array[currentPoint[0] + o[0], currentPoint[1] + o[1]].temperature
                maxTempSolid = max(t, maxTempSolid) if self.array[currentPoint[0] + o[0], currentPoint[1] + o[1]].material in MaterialType.SOLID else maxTempSolid

            r = self.array[currentPoint].r
            area: pint.Quantity = self.array[currentPoint].area

            score = (maxTempSolid.m_as("degR") / (area.m_as("in^2")**2) * r**2)**1.5
            pointData.append((currentPoint, score))
            totalScore += score

        # should now have a total score and each point contributes to that
        # the deltaT of each point will be directly proportional to its scores percentage

        runningDt = self.array[startPoint].temperature
        for point, score in pointData:
            runningDt += score/totalScore * deltaT
            self.array[point].temperature = runningDt

        longLength = 0
        shortLength = 0

        for point in alive_it(channelPoints[::-1]): # start pressure at outlet
            if self.array[point].material == DomainMaterial.COOLANT_OUTLET:
                continue

            futureFlow = self.array[point].futureFlow
            coolingLoopData = self.coolingLoops[self.array[point].id]
            futurePress = self.array[futureFlow].pressure
            deltaL = Q_(np.sqrt((self.array[point].x - self.array[futureFlow].x)**2 + (self.array[point].r - self.array[futureFlow].r)**2), unitReg.inch).to(unitReg.ft)
            (mu, cp, _, _, rho, _, _, _, _) = fluid.get_fluid_properties(coolingLoopData.fluid, self.array[futureFlow].temperature, self.array[futureFlow].pressure)
            mdotperchannel = coolingLoopData.mdot / coolingLoopData.numChannels
            Re = mdotperchannel*self.array[futureFlow].hydraulicDiameter/mu/self.array[futureFlow].area   # Reynolds number
            Re = Re.to(unitReg.dimensionless)
            f = fluid.DarcyFrictionFactor(Re, Q_(.05, unitReg.inch), self.array[futureFlow].hydraulicDiameter)
            vel = (mdotperchannel / rho / self.array[futureFlow].area).to(unitReg.ft/unitReg.sec)
            dp = fluid.FrictionPressureLoss(f, deltaL, self.array[futureFlow].hydraulicDiameter, rho, vel).to(unitReg.psi)
            self.array[point].pressure = futurePress + dp
            self.array[point].velocity = vel.to(unitReg.ft/unitReg.sec)

            if point[0] != 27:
                shortLength += deltaL
            else:
                longLength += deltaL

        for i in range(self.vpoints):
            for j in range(self.hpoints):
                if self.array[i,j].material == DomainMaterial.COOLANT_BULK and self.array[i,j].id == loopID:
                    self.array[i,j].temperature = self.array[self.array[i,j].previousFlow].temperature
                    self.array[i,j].pressure = self.array[self.array[i,j].previousFlow].pressure
                    self.array[i,j].velocity = self.array[self.array[i,j].previousFlow].velocity

        # print(f"long length: {longLength.to(unitReg.inch)}")
        # print(f"short length: {shortLength.to(unitReg.inch)}")


        print("done????")

    def AssignChamberTemps(self, chamber: np.ndarray, exhaust: Gas, startPoint: tuple, endPoint: tuple, chamberWallRadius: pint.Quantity, plugBase: pint.Quantity, Astar: pint.Quantity):
        tic = time.perf_counter()
        print("assigning stagnant")
        flowAngle = np.arctan((endPoint[1] - startPoint[1])/(endPoint[0] - startPoint[0]))
        phi = np.pi/2 - flowAngle

        farAway1 = (startPoint[0] - chamberWallRadius.magnitude*np.sin(flowAngle), startPoint[1] + chamberWallRadius.magnitude*np.cos(flowAngle))
        farAway2 = (startPoint[0] + chamberWallRadius.magnitude*np.sin(flowAngle), startPoint[1] - chamberWallRadius.magnitude*np.cos(flowAngle))
        startPointU, _ = material.intersectPolyAt(chamber, startPoint, farAway1)
        startPointL, _ = material.intersectPolyAt(chamber, startPoint, farAway2)

        originX = np.linspace(startPoint[0], endPoint[0], int(self.hpoints))
        originR = np.linspace(startPoint[1], endPoint[1], int(self.hpoints))

        iStart, jStart = self.ChamberStartCell()

        for j in range(jStart, self.hpoints):
            for i in range(iStart, self.vpoints):
                if self.array[i,j].material != DomainMaterial.CHAMBER:
                    break
                if self.lineInCell(startPointU, startPointL, i, j):
                    jStart = j if i == iStart else jStart
                    continue
                self.array[i,j].temperature = exhaust.stagTemp.to(unitReg.degR)
                self.array[i,j].pressure = exhaust.stagPress.to(unitReg.psi)
                self.array[i,j].velocity = Q_(1, unitReg.foot/unitReg.sec)
                self.array[i,j].hydraulicDiameter = 2*(chamberWallRadius - plugBase).to(unitReg.inch)
                self.array[i,j].area = Q_(np.pi*(chamberWallRadius - plugBase)**2, unitReg.inch**2)
                self.array[i,j].flowHeight = Q_(chamberWallRadius - plugBase, unitReg.inch)
            if self.lineInCell(startPointU, startPointL, i, j):
                break

        print(f"Assinging straight flow")
        for oX, oR in alive_it(zip(originX, originR)):
            farAway1 = (oX - chamberWallRadius.magnitude*np.sin(flowAngle), oR + chamberWallRadius.magnitude*np.cos(flowAngle))
            farAway2 = (oX + chamberWallRadius.magnitude*np.sin(flowAngle), oR - chamberWallRadius.magnitude*np.cos(flowAngle))
            startPointU, _ = material.intersectPolyAt(chamber, (oX, oR), farAway1)
            startPointL, _ = material.intersectPolyAt(chamber, (oX, oR), farAway2)

            # plt.plot([startPointU[0], startPointL[0]], [startPointU[1], startPointL[1]], '-gx')

            area = np.pi/np.sin(phi) * (startPointU[1]**2 - startPointL[1]**2)
            height = np.sqrt((startPointU[0] - startPointL[0])**2 + (startPointU[1] - startPointL[1])**2)

            AAstar = Q_(area, unitReg.inch**2)/Astar
            mach = fsolve(lambda M: gas.Isentropic1DExpansion(M, exhaust.gammaTyp) - AAstar, .25)[0]
            temperature = (gas.StagTempRatio(mach, exhaust) * exhaust.stagTemp)
            velocity = mach * np.sqrt(exhaust.getVariableGamma(mach) * exhaust.Rgas * temperature)
            pressure = gas.StagPressRatio(mach, exhaust) * exhaust.stagPress
            hydroD = Q_(2*np.sqrt((startPointU[0] - startPointL[0])**2 + (startPointU[1] - startPointL[1])**2), unitReg.inch)

            cells = self.cellsOnLine(startPointL, startPointU)
            for i,j in cells:
                if self.array[i,j].material == DomainMaterial.CHAMBER:
                    self.array[i,j].temperature = temperature.to(unitReg.degR)
                    self.array[i,j].velocity = velocity.to(unitReg.ft/unitReg.s)
                    self.array[i,j].pressure = pressure.to(unitReg.psi)
                    self.array[i,j].hydraulicDiameter = hydroD.to(unitReg.inch)
                    self.array[i,j].area = Q_(area, unitReg.inch**2)
                    self.array[i,j].flowHeight = Q_(height, unitReg.inch)


        # curve section
        print(f"assigning bend flow")
        farAway = (endPoint[0] - chamberWallRadius.magnitude*np.sin(flowAngle), endPoint[1] - chamberWallRadius.magnitude*np.cos(flowAngle))
        plt.plot(farAway[0], farAway[1], 'rx')
        
        intersect, indexPair = material.intersectPolyAt(chamber, endPoint, farAway)
        distance = np.sqrt((chamber[indexPair[0]].x - chamber[indexPair[1]].x)**2 + (chamber[indexPair[0]].r - chamber[indexPair[1]].r)**2)
        interpFactor = int(np.ceil(5*distance / self.xstep))
        ic(interpFactor)
        plt.plot([intersect[0]], [intersect[1]], 'gx')

        iiContourStart = indexPair[0]

        curAngle = np.pi/2 + np.arctan2(chamber[iiContourStart+2].r - chamber[iiContourStart+1].r, chamber[iiContourStart+2].x - chamber[iiContourStart+1].x)
        prevAngle = np.pi/2 + np.arctan2(chamber[iiContourStart+1].r - chamber[iiContourStart].r, chamber[iiContourStart+1].x - chamber[iiContourStart].x)
        angle0 = prevAngle
        ii = iiContourStart + 1
        ii1 = ii + 1
        # cArr = ['.b-', '.k-', '.k-', '.k-', '.k-', '.k-', '.k-', '.k-', '.k-', '.k-', '.k-', '.k-']
        # ic(angle)
        while curAngle <= angle0 + 1e-3:
            # ci = 0
            angles = np.linspace(prevAngle, curAngle, interpFactor + 1)
            points = np.linspace((chamber[ii].x, chamber[ii].r), (chamber[ii+1].x, chamber[ii+1].r), interpFactor + 1)
            for angle, lowerPoint in zip(angles[:-1], points[:-1]):
                farAway = (lowerPoint[0] + chamberWallRadius.magnitude*np.cos(angle), lowerPoint[1] + chamberWallRadius.magnitude*np.sin(angle))
                upperPoint, _ = material.intersectPolyAt(chamber, (lowerPoint[0], lowerPoint[1] + 1e-3), farAway)

                # plt.plot([lowerPoint[0], upperPoint[0]], [lowerPoint[1], upperPoint[1]], cArr[ci % 10], linewidth=.5)
                # ci += 1

                area = np.pi/np.sin(angle) * (upperPoint[1]**2 - lowerPoint[1]**2)
                height = np.sqrt((upperPoint[0] - lowerPoint[0])**2 + (upperPoint[1] - lowerPoint[1])**2)

                AAstar = Q_(area, unitReg.inch**2)/Astar
                if AAstar < 1:
                    mach = 1
                else:
                    mach = fsolve(lambda M: gas.Isentropic1DExpansion(M, exhaust.gammaTyp) - AAstar, .25)[0]
                temperature = (gas.StagTempRatio(mach, exhaust) * exhaust.stagTemp)
                pressure = gas.StagPressRatio(mach, exhaust) * exhaust.stagPress
                velocity = mach * np.sqrt(exhaust.getVariableGamma(mach) * exhaust.Rgas * temperature)
                hydroD = Q_(2*np.sqrt((upperPoint[0] - lowerPoint[0])**2 + (upperPoint[1] - lowerPoint[1])**2), unitReg.inch)

                cells = self.cellsOnLine(lowerPoint, upperPoint)
                for i,j in cells:
                    if self.array[i,j].material == DomainMaterial.CHAMBER:
                        self.array[i,j].temperature = temperature.to(unitReg.degR)
                        self.array[i,j].velocity = velocity.to(unitReg.ft/unitReg.s)
                        self.array[i,j].pressure = pressure.to(unitReg.psi)
                        self.array[i,j].hydraulicDiameter = hydroD.to(unitReg.inch)
                        self.array[i,j].area = Q_(area, unitReg.inch**2)
                        self.array[i,j].flowHeight = Q_(height, unitReg.inch)

            ii = ii1

            ii1 = ii + 1
            if ii1 >= len(chamber):
                ii1 = 0
            prevAngle = curAngle
            curAngle = np.pi/2 + np.arctan2(chamber[ii1].r - chamber[ii].r, chamber[ii1].x - chamber[ii].x)
        
        toc = time.perf_counter()
        print(f"Time to assign chamber temps: {toc - tic}")

    def AssignExternalTemps(self, gridField: np.ndarray[CharacteristicPoint], exhaust: Gas, Astar: pint.Quantity):
        points = []
        values = []
        minc = (2e15, 2e15)
        maxc = (-2e15, -2e15)
        for row in gridField:
            for point in row:
                if point is None:
                    continue
                if point.x < minc[0]:
                    minc = (point.x, minc[1])
                if point.x > maxc[0]:
                    maxc = (point.x, maxc[1])
                if point.r < minc[1]:
                    minc = (minc[0], point.r)
                if point.r > maxc[1]:
                    maxc = (maxc[0], point.r)
                points.append((point.x, point.r))
                if isinstance(point, pint.Quantity):
                    values.append(point.mach.magnitude)
                else:
                    values.append(point.mach)

        minc = (max(minc[0], self.x0), max(minc[1], self.r0 - self.height))
        maxc = (min(maxc[0], self.x0 + self.width), min(maxc[1], self.r0))
        mincell = self.CoordsToCell(minc[0], minc[1])
        maxcell = self.CoordsToCell(maxc[0], maxc[1])
        rows, cols = np.mgrid[maxcell[0]:mincell[0], mincell[1]:maxcell[1]]
        x = np.array([[p.x for p in row] for row in self.array[rows, cols]])
        r = np.array([[p.r for p in row] for row in self.array[rows, cols]])

        machs = griddata(points, values, (x,r))

        for i in range(np.size(machs, 0)):
            for j in range(np.size(machs, 1)):
                if np.isnan(machs[i,j]):
                    continue
                mach = machs[i,j]
                rowDomain = rows[i,j]
                colDomain = cols[i,j]
                self.array[rowDomain, colDomain].material = DomainMaterial.CHAMBER
                self.array[rowDomain, colDomain].temperature = gas.StagTempRatio(mach, exhaust) * exhaust.stagTemp
                self.array[rowDomain, colDomain].velocity = mach * np.sqrt(exhaust.getVariableGamma(mach) * exhaust.Rgas * self.array[rowDomain, colDomain].temperature)
                self.array[rowDomain, colDomain].area = gas.Isentropic1DExpansion(mach, exhaust.gammaTyp) * Astar
                

        
        print("done")
        
    def AssignBorders(self, fluidID = -1):
        finalI = self.vpoints - 1
        finalJ = self.hpoints - 1

        with alive_bar(self.vpoints*self.hpoints) as bar:
            for i in range(self.vpoints):
                for j in range(self.hpoints):
                    if fluidID != -1 and self.array[i,j].id != fluidID:
                        bar()
                        continue
                    lowI = max(i - 1, 0)
                    lowJ = max(j - 1, 0)
                    highI = min(i + 1, finalI)
                    highJ = min(j + 1, finalJ)
                    checks = np.array([self.array[lowI, j].material, self.array[highI, j].material, self.array[i, lowJ].material, self.array[i, highJ].material])
                    equals = checks != self.array[i,j].material
                    coolantWalls = [check in MaterialType.COOLANT_WALL for check in checks]
                    self.array[i, j].border = np.any(equals ^ coolantWalls) if self.array[i,j].material in MaterialType.COOLANT else np.any(equals)
                    bar()

    def cellsOnLine(self, point1, point2):
        linedx = point2[0] - point1[0]
        linedr = point2[1] - point1[1]
        cellsx = int(linedx / self.xstep)
        cellsr = int(linedr / self.rstep)

        cellsLine = max(int(np.sqrt(cellsx**2 + cellsr**2) * 1.5),10)
        x = np.linspace(point1[0], point2[0], cellsLine)
        r = np.linspace(point1[1], point2[1], cellsLine)

        return [self.CoordsToCell(xp, rp) for xp, rp in zip(x,r)]    

    def isInCell(self, point, row, col):
        xmin = self.array[row, col].x - self.xstep/2
        xmax = self.array[row, col].x + self.xstep/2
        ymin = self.array[row, col].r - self.rstep/2
        ymax = self.array[row, col].r + self.rstep/2
        return xmin <= point[0] and point[0] <= xmax and ymin <= point[1] and point[1] <= ymax
        
    def lineInCell(self, point1, point2, row, col, res = -1):
        m = (point2[1] - point1[1]) / (point2[0] - point1[0])
        a = -m
        b = 1
        c = m*point1[0] - point1[1]
        x0 = self.array[row, col].x
        y0 = self.array[row, col].r

        x = (b*(b*x0 - a*y0) - a*c)/(a**2 + b**2)
        y = (a*(-b*x0 + a*y0) - b*c)/(a**2 + b**2)
        return self.isInCell((x,y), row, col)
    
    def ChamberStartCell(self):
        for i in range(self.vpoints):
            for j in range(self.hpoints):
                if self.array[i,j].material == DomainMaterial.CHAMBER:
                    return (i, j)

    def CoordsToCell(self, x, r):
        dx = x - self.x0 + self.xstep/2
        dr = self.r0 - r + self.rstep/2
        i = int((dr/self.rstep))
        j = int((dx/self.xstep))
        return (i, j)

    def DumpFile(self, filename):
        joblib.dump(self, filename + '.msh.z', compress=True)

    def ApplyStateMap(self, prevDomain: 'DomainMC', states: set):
        for i in range(self.vpoints):
            for j in range(self.hpoints):
                if self.array[i,j].material in material.MaterialType.STATIC_TEMP:
                    continue
                refPoint = prevDomain.CoordsToCell(self.array[i,j].x, self.array[i,j].r)
                if prevDomain.array[refPoint].material != self.array[i,j].material:
                    ir = refPoint[0]
                    jr = refPoint[1]
                    offsets = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
                    minDist = 1e15
                    for o in offsets:
                        if 0 <= ir + o[0] < prevDomain.vpoints and 0 <= jr + o[1] < prevDomain.hpoints:
                            d = np.sqrt((self.array[i,j].x - prevDomain.array[ir + o[0], jr + o[1]].x)**2 + (self.array[i,j].r - prevDomain.array[ir + o[0], jr + o[1]].r)**2)
                            if prevDomain.array[ir + o[0], jr + o[1]].material == self.array[i,j].material and d < minDist:
                                refPoint = (ir + o[0], jr + o[1])
                                minDist = d
                for state in states:
                    self.array[i,j].__setattr__(state, prevDomain.array[refPoint].__getattribute__(state))

    @staticmethod
    def LoadFile(filename):
        print("Loading file")
        loaded: DomainMC = joblib.load(filename + '.msh.z')
        # DomainMC.ConvertUnits(loaded)
        return loaded
    
    @staticmethod
    def ConvertUnits(domain):
        print("Converting units")
        attributes = list(domain.array[0,0].__dict__.keys())
        for attr in attributes:
            testAttr = domain.array[0,0].__getattribute__(attr)
            if isinstance(testAttr, pint.Quantity):
                with alive_bar(domain.vpoints*domain.hpoints) as bar:
                    print(f"Converting {attr}")
                    unit = str(testAttr.units)
                    t = [[domain.array[i,j].__getattribute__(attr).magnitude for i in range(domain.vpoints)] for j in range(domain.hpoints)]
                    for i in range(domain.vpoints):
                        for j in range(domain.hpoints):
                            domain.array[i,j].__setattr__(attr, Q_(t[j][i], unit))
                            bar()
        return domain

def EvalMaterialProcess2(i, hsteps, pn, gridData, size, cowl, chamber, plug):
    res = []

    x0 = gridData[0]
    xstep = gridData[1]
    r0 = gridData[2]
    rstep = gridData[3]

    for j in range(hsteps):
        point = (x0 + j*xstep, r0 - i*rstep)
        material = AssignMaterial(point, size, np.array([]), cowl, chamber, plug)
        res.append((i, j, material))
    return res

def AssignMaterial(point, size, coolant, cowl, chamber, plug):
    if material.isIntersect(point, coolant, size):
        return DomainMaterial.COOLANT
    if material.isIntersect(point, cowl, size):
        return DomainMaterial.COWL
    if material.isIntersect(point, chamber, size):
        return DomainMaterial.CHAMBER
    if material.isIntersect(point, plug, size):
        return DomainMaterial.PLUG
    return DomainMaterial.FREE