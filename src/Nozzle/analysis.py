import numpy as np
np.product = np.prod
from dataclasses import dataclass
import matplotlib.pyplot as plt
import fluids.gas as gas
from fluids.gas import MachAngle, mach2machStar, machStar2mach, Gas
import General.design as DESIGN
import Nozzle.nozzle as nozzle
from icecream import ic
import Nozzle.config as config
from General.units import Q_, unitReg
import matrix_viewer as mv
import logging
from General.setenv import setupLogging
setupLogging()

@dataclass
class CharacteristicPoint:
    x: float
    r: float
    theta: float
    machStar: float
    s: float = 0
    mach: float = 0
    alpha: float = 0
    terminate: bool = False

    F: float = 0
    G: float = 0
    H: float = 0
    J: float = 0

    def __repr__(self) -> str:
        return f"CP @ ({self.x:.3f}, {self.r:.3f}), theta: {np.rad2deg(self.theta):.3f}, s: {self.s:.3f}, mach: {self.mach:.3f}"

    def clone(self):
        return CharacteristicPoint(self.x, self.r, self.theta, self.machStar, self.s, self.mach, self.alpha, terminate=self.terminate)

    def setTerminate(self):
        self.terminate = True
        return self

    def CalculateRightVariant(self, gas: Gas) -> 'CharacteristicPoint': # I characteristic
        Rgas = gas.Rgas.magnitude
        gamma = gas.gammaTyp
        RV = self.clone()
        RV.F = np.tan(RV.theta - RV.alpha)
        RV.G = 1/np.tan(RV.alpha)/RV.machStar
        RV.H = -np.sin(RV.theta)*np.sin(RV.alpha)/(RV.r*np.sin(RV.theta - RV.alpha))
        RV.J = np.sin(RV.alpha)*np.cos(RV.alpha)/(Rgas * gamma)
        return RV

    def CalculateLeftVariant(self, gas: Gas) -> 'CharacteristicPoint': # II characteristic
        Rgas = gas.Rgas.magnitude
        gamma = gas.gammaTyp
        LV = self.clone()
        LV.F = np.tan(LV.theta + LV.alpha)
        LV.G = -1/np.tan(LV.alpha)/LV.machStar
        LV.H = np.sin(LV.theta)*np.sin(LV.alpha)/(LV.r*np.cos(LV.theta + LV.alpha)) # makes L.H nan from LV.r
        LV.J = -np.sin(LV.alpha)*np.cos(LV.alpha)/(Rgas * gamma)
        return LV
    
    def setCoefficients(self, F, G, H, J) -> None:
        self.F = F
        self.G = G
        self.H = H
        self.J = J

    def NextIterationPoint(self, N: 'CharacteristicPoint') -> 'CharacteristicPoint':
        newF = (self.F + N.F)/2
        newG = (self.G + N.G)/2
        newH = (self.H + N.H)/2
        newJ = (self.J + N.J)/2
        nextPoint = self.clone()
        nextPoint.setCoefficients(newF, newG, newH, newJ)

        return nextPoint

    # Field Point Calculations
    @staticmethod
    def CalculateSolidBoundaryIntersect(point: 'CharacteristicPoint', contour: np.ndarray[nozzle.ContourPoint], isRight: bool) -> tuple[float, float] | None:        
        Sx, Sy = point.x, point.r
        angle = point.theta - point.alpha if isRight else point.theta + point.alpha

        for j in range(len(contour) - 1):
            a, b, c, d = contour[j].x, contour[j].r, contour[j+1].x, contour[j+1].r
            mL = np.tan(angle)
            mC = (d - b)/(c - a)
            Amat = np.array([[mL, -1], [mC, -1]])
            bmat = np.array([[mL*Sx - Sy], [mC*a - b]])
            X = np.linalg.solve(Amat, bmat)
            if a <= X[0,0] <= c:
                Bx = X[0,0]
                By = X[1,0]
                # plt.plot([a, c], [b, d], '-or', linewidth=1)
                return (Bx, By), np.atan2((d - b),(c - a))
        return None, None

    @staticmethod
    def CalculateFieldPoint(L: 'CharacteristicPoint', R: 'CharacteristicPoint', workingGas: Gas, tol = 1e-6) -> 'CharacteristicPoint':
        L = L.CalculateLeftVariant(workingGas)
        R = R.CalculateRightVariant(workingGas)
        N = CharacteristicPoint.ApproxCharacteristicEqn(L, R, workingGas)

        for i in range(30):
            NR = N.CalculateRightVariant(workingGas)
            NL = N.CalculateLeftVariant(workingGas)

            L2 = L.NextIterationPoint(NL)
            R2 = R.NextIterationPoint(NR)

            NN = CharacteristicPoint.ApproxCharacteristicEqn(L2, R2, workingGas)
            if np.abs((NN.theta - N.theta)/(NN.theta)) < tol:
                logging.debug(f"Converged in {i} iterations")
                return NN
            else:
                N = NN.clone()
                del NN, NL, NR, L2, R2

        logging.debug(f"Did not converge in 30 iterations")
        return N

    @staticmethod
    def ApproxCharacteristicEqn(L: 'CharacteristicPoint', R: 'CharacteristicPoint', workingGas: Gas) -> 'CharacteristicPoint':
        Amat = np.array([[1, -R.F], [1, -L.F]])
        b = np.array([[R.r - R.F*R.x], [L.r - L.F*L.x]])
        X = np.linalg.solve(Amat, b)
        r = X[0, 0]
        x = X[1, 0]

        nr = (x - R.x)*np.sin(R.alpha)/np.cos(R.theta - R.alpha)
        nl = (x - L.x)*np.sin(L.alpha)/np.cos(L.theta + L.alpha)
        s = (R.s - L.s) / (nl + nr)
        s = R.s + s*nr

        if abs(L.r) > 1e-2: # L is not on axis
            Amat = np.array([[1, R.G], [1, L.G]])
            b = np.array([[R.theta + R.G*R.machStar - R.H*(r - R.r) - R.J*(s - R.s)], [L.theta + L.G*L.machStar - L.H*(x - L.x) - L.J*(s - L.s)]])
        else:
            Amat = np.array([[1, R.G], [2, L.G]])
            b = np.array([[R.theta + R.G*R.machStar - R.H*(r - R.r) - R.J*(s - R.s)], [L.theta + L.G*L.machStar - L.J*(s - L.s)]])
        
        X = np.linalg.solve(Amat, b)
        theta = X[0, 0]
        machStar = X[1, 0]
        mach = machStar2mach(machStar, workingGas.gammaTyp)
        alpha = MachAngle(mach)

        return CharacteristicPoint(x, r, theta, machStar, s, mach, alpha)

    @staticmethod
    def CalculateSolidReflect(point: 'CharacteristicPoint', isRight: bool, contour: np.ndarray[nozzle.ContourPoint], workingGas: Gas, tol = 1e-6) -> 'CharacteristicPoint':
        PV = point.CalculateRightVariant(workingGas) if isRight else point.CalculateLeftVariant(workingGas)
        N = CharacteristicPoint.ApproxSolidReflect(PV, isRight, contour, workingGas)
        if N is None:
            return None
        for i in range(30):
            N2 = N.CalculateRightVariant(workingGas) if isRight else N.CalculateLeftVariant(workingGas)
            N2 = PV.NextIterationPoint(N2)

            NN = CharacteristicPoint.ApproxSolidReflect(N2, isRight, contour, workingGas)
            if NN is None:
                return None
            if np.abs((NN.theta - N.theta)/(NN.theta)) < tol:
                logging.debug(f"Solid Converged in {i} iterations")
                return NN
            else:
                N = NN.clone()
                del NN, N2

        logging.debug(f"Solid Did not converge in 30 iterations")
        return N

    @staticmethod
    def ApproxSolidReflect(PV: 'CharacteristicPoint', isRight: bool, contour: np.ndarray[nozzle.ContourPoint], workingGas: Gas) -> 'CharacteristicPoint':
        # ic(PV)
        intersect, theta = CharacteristicPoint.CalculateSolidBoundaryIntersect(PV, contour, isRight)
        if intersect is None:
            return None#CharacteristicPoint.CalculateGasReflect(point, isRight, PbPc, workingGas)
        
        x = intersect[0]
        r = intersect[1]
        s = 0 # TODO entropy of contour streamline

        Hfact = (r - PV.r) if isRight else (x - PV.x)

        machStar = PV.machStar + (-(theta - PV.theta) - PV.H*Hfact - PV.J*(s - PV.s))/PV.G

        return CharacteristicPoint(x, r, theta, machStar, s, machStar2mach(machStar, workingGas.gammaTyp), MachAngle(machStar2mach(machStar, workingGas.gammaTyp)))

    @staticmethod
    def CalculateGasReflect(point: 'CharacteristicPoint', isRight: bool, PambPc: float, workingGas: Gas, streamline: np.ndarray['CharacteristicPoint'], tol = 1e-6) -> 'CharacteristicPoint':
        PV = point.CalculateRightVariant(workingGas) if isRight else point.CalculateLeftVariant(workingGas)
        N = CharacteristicPoint.ApproxGasReflect(PV, isRight, PambPc, workingGas, streamline)

        for i in range(30):
            N2 = N.CalculateRightVariant(workingGas) if isRight else N.CalculateLeftVariant(workingGas)
            N2 = PV.NextIterationPoint(N2)

            NN = CharacteristicPoint.ApproxGasReflect(N2, isRight, PambPc, workingGas, streamline)
            if np.abs((NN.theta - N.theta)/(NN.theta)) < tol:
                logging.debug(f"Gas Converged in {i} iterations")
                return NN
            else:
                N = NN.clone()
                del NN, N2

        logging.debug(f"Gas Did not converge in 30 iterations")
        return N

    @staticmethod
    def ApproxGasReflect(PV: 'CharacteristicPoint', isRight, PambPc, workingGas: Gas, streamline: np.ndarray['CharacteristicPoint']) -> 'CharacteristicPoint':
        machInf = np.sqrt((PambPc**(-1/workingGas.gammaTyp[5]) - 1)/workingGas.gammaTyp[2])
        streamPoint: CharacteristicPoint = streamline[-1]

        angle = PV.theta - PV.alpha if isRight else PV.theta + PV.alpha
        mC = np.tan(angle)
        mS = np.tan(streamPoint.theta)
        Amat = np.array([[mC, -1], [mS, -1]])
        bmat = np.array([[mC*PV.x - PV.r], [mS*streamPoint.x - streamPoint.r]])
        X = np.linalg.solve(Amat, bmat)
        x = X[0, 0]
        r = X[1, 0]
        s = streamPoint.s
        machStar = mach2machStar(machInf, workingGas.gammaTyp)
        Hfact = (r - PV.r) if isRight else (x - PV.x)
        theta = PV.theta - PV.G*(machStar - PV.machStar) - PV.H*Hfact - PV.J*(s - PV.s)

        return CharacteristicPoint(x, r, theta, machStar, s, machStar2mach(machStar, workingGas.gammaTyp), MachAngle(machStar2mach(machStar, workingGas.gammaTyp)))
    
    @staticmethod
    def CaclulateAxisReflect(point: 'CharacteristicPoint', isRight: bool, workingGas, tol = 1e-6) -> 'CharacteristicPoint':
        PV = point.CalculateRightVariant(workingGas) if isRight else point.CalculateLeftVariant(workingGas)
        N = CharacteristicPoint.ApproxAxisReflect(PV, isRight, workingGas)

        return N

        for i in range(30):
            N2 = N.CalculateRightVariant(workingGas) if isRight else N.CalculateLeftVariant(workingGas)
            N2 = PV.NextIterationPoint(N2)

            NN = CharacteristicPoint.ApproxAxisReflect(N2, isRight, workingGas)
            if np.abs((NN.mach - N.mach)/(NN.mach)) < tol: # mach not theta because theta is 0
                logging.debug(f"Axis converged in {i} iterations")
                return NN
            else:
                N = NN.clone()
                del NN, N2

        logging.debug(f"Axis did not converge in 30 iterations")
        return N
    
    @staticmethod
    def ApproxAxisReflect(PV: 'CharacteristicPoint', isRight: bool, workingGas: Gas) -> 'CharacteristicPoint':
        r = 0
        theta = 0

        x = (r - PV.r)/PV.F + PV.x

        s = 0 # TODO entropy of contour streamline
        Hfact = (r - PV.r) if isRight else (x - PV.x)
        machStar = (-(theta - PV.theta) - PV.H*Hfact - PV.J*(s - PV.s))/PV.G + PV.machStar

        return CharacteristicPoint(x, r, theta, machStar, s, machStar2mach(machStar, workingGas.gammaTyp), MachAngle(machStar2mach(machStar, workingGas.gammaTyp)))

def CalculateComplexField(contour, Pamb: Q_, workingGas: Gas, Mt: float, Tt: float, Rt: Q_, scale = 1, Rsteps = 20, Lsteps = 0, reflections = 3, fig = None):
    PbPc = DESIGN.basePressure/DESIGN.chamberPressure
    PambPc = Pamb/DESIGN.chamberPressure
    gamma = workingGas.gammaTyp
    Rt = Rt.to(unitReg.inch).magnitude
    xt = (scale - Rt)*np.tan(Tt)
    Me = np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2])
    thetaExit = Tt + gas.PrandtlMeyerFunction(Me, gamma) - gas.PrandtlMeyerFunction(Mt, gamma)

    outerStreamLine = np.array([CharacteristicPoint(0, scale, thetaExit, mach2machStar(Me, gamma), 0, Me, MachAngle(Me))])
    innerStreamLine = np.array([CharacteristicPoint(xt, Rt, Tt, mach2machStar(Mt, gamma), 0, Mt, MachAngle(Mt))])

    rLines = np.empty((Rsteps, 1 + (Lsteps + Rsteps)*reflections), dtype=CharacteristicPoint)
    lLines = np.empty((Lsteps, 1 + (Lsteps + Rsteps)*reflections), dtype=CharacteristicPoint)

    rLines[:, 0] = np.transpose(GenerateExpansionFan(Me, Mt, Tt, workingGas, Rsteps, scale))
    lLines[:, 0] = np.transpose(GenerateStartLine(Rt, Mt, Tt, workingGas, Lsteps, scale))

    streamlines = (innerStreamLine, outerStreamLine)

    # plt.ion()
    # plt.show()
    for i in range(reflections):
        rLines, lLines = PropogateRegionAll(rLines, lLines, workingGas, i)
        rLines, lLines, streamlines = ReflectionRegionAll(rLines, lLines, contour, PambPc, PbPc, streamlines, workingGas, i, fig)

    return rLines, lLines, streamlines

def GenerateStartLine(Rt: float, machT, thetaT, workingGas: Gas, arraySize: int, scale = 1):
    xt: Q_ = (scale - Rt)*np.tan(thetaT)

    x = np.linspace(xt, 0, arraySize + 2)
    r = np.linspace(Rt, scale, arraySize + 2)
    startline = np.array([CharacteristicPoint(x[i], r[i], thetaT, mach2machStar(machT, workingGas.gammaTyp), 0, machT, MachAngle(machT)) for i in range(1, arraySize+1)])
    return startline[::-1]

def GenerateExpansionFan(machE: float, machT: float, thetaT: float, workingGas: Gas, arraySize: int, scale = 1):
    gamma = workingGas.gammaTyp
    machs = np.linspace(machT, machE, arraySize) if machT > config.MIN_MOC_MACH else np.linspace(config.MIN_MOC_MACH, machE, arraySize)

    thetas = thetaT + gas.PrandtlMeyerFunction(machs, gamma) - gas.PrandtlMeyerFunction(machT, gamma)

    expansionFanArray = np.empty(arraySize, dtype=CharacteristicPoint)
    for i, mach in enumerate(machs):
        expansionFanArray[i] = CharacteristicPoint(0, scale, thetas[i], mach2machStar(mach, gamma), 0, mach, MachAngle(mach))

    return expansionFanArray

def PropogateRegionAll(rLines: np.ndarray[CharacteristicPoint], lLines: np.ndarray[CharacteristicPoint], workingGas: Gas, reflection: int):
    R0: int = rLines.shape[0]
    L0: int = lLines.shape[0]

    start = 1 + (reflection)*(R0 + L0) # reflection should start at 0
    region = np.empty((R0+1, L0+1), dtype=CharacteristicPoint)
    region[1:,0] = rLines[:,start-1]
    region[0,1:] = np.transpose(lLines[:,start-1])
    for i in range(1, R0+1):
        for j in range(1, L0+1):
            region[i,j] = CharacteristicPoint.CalculateFieldPoint(region[i-1,j], region[i,j-1], workingGas) if reflection % 2 == 0 else CharacteristicPoint.CalculateFieldPoint(region[i,j-1], region[i-1,j], workingGas)
    
    rLines[:,start:start+L0] = region[1:,1:]
    lLines[:,start:start+R0] = np.transpose(region[1:,1:])

    return rLines, lLines

def ReflectionRegionAll(rLines: np.ndarray[CharacteristicPoint], lLines: np.ndarray[CharacteristicPoint], contour, PambPc, PbPc, streamlines, workingGas: Gas, reflection: int, fig):
    R0: int = rLines.shape[0]
    L0: int = lLines.shape[0]
    innerStreamline, outerStreamline = streamlines
    
    rlines, innerStreamline, outerStreamline = ReflectionRegion(rLines, R0, L0, contour, PambPc, PbPc, innerStreamline, outerStreamline, workingGas, reflection, True, fig)
    llines, innerStreamline, outerStreamline = ReflectionRegion(lLines, R0, L0, contour, PambPc, PbPc, innerStreamline, outerStreamline, workingGas, reflection, False, fig)

    return rlines, llines, (innerStreamline, outerStreamline)

def ReflectionRegion(lines: np.ndarray[CharacteristicPoint], R0, L0, contour, PambPc, PbPc, innerStreamline, outerStreamline, workingGas: Gas, reflection, startAsRight: bool, fig): # startAsRight is true if region is being calculated in rlines
    X0 = R0 if startAsRight else L0
    Y0 = L0 if startAsRight else R0
    start = 1 + Y0 + (reflection)*(X0 + Y0) # reflection should start at 0
    region = np.empty((X0+1, X0+1), dtype=CharacteristicPoint)
    region[1:,0] = lines[:,start-1]
    region[0,1:] = np.transpose(lines[:,start-1])

    for i in range(1, X0+1):
        for j in range(1, i+1):
            # ic(i,j)
            # PlotCharacteristicLines(fig, region[1:, 1:])
            # fig.canvas.draw()
            # fig.canvas.flush_events()
            # plt.waitforbuttonpress()
            # fig.axes[0].clear()
            isRight = not (reflection % 2 == 0) ^ startAsRight
            if region[i, j-1].terminate:
                region[i,j] = region[i-1,j].clone()
                continue
            if i == j:
                reflectOrigin = region[i, j-1] # previous point in the same line
                region[i,j], innerStreamline, outerStreamline = DoReflect(reflectOrigin, isRight, innerStreamline, outerStreamline, contour, PambPc, PbPc, workingGas)
            else:
                region[i,j] = CharacteristicPoint.CalculateFieldPoint(region[i-1,j], region[i,j-1], workingGas) if isRight else CharacteristicPoint.CalculateFieldPoint(region[i,j-1], region[i-1,j], workingGas)
                region[j, i] = region[i, j].clone()

    lines[:, start:start+X0] = region[1:,1:]
    # plt.waitforbuttonpress()
    return lines, innerStreamline, outerStreamline

def DoReflect(point: CharacteristicPoint, isRight: bool, innerStreamline: np.ndarray, outerStreamline: np.ndarray, contour: np.ndarray, PambPc: float, PbPc: float, workingGas: Gas):
    if isRight:
        newPoint = None
        doInnerLine = innerStreamline[-1].r > 1e-3
        doIntersect = innerStreamline[-1].r > contour[-1].r
        # ic(innerStreamline[-1].r, doInnerLine)
        if doInnerLine and doIntersect:
            newPoint = CharacteristicPoint.CalculateSolidReflect(point, isRight, contour, workingGas)
        if newPoint is None and doInnerLine: # it missed the contour
            if doIntersect:
                newAngle = np.atan2(contour[-1].r - innerStreamline[-1].r, contour[-1].x - innerStreamline[-1].x)
                innerStreamline[-1].theta = newAngle
            newPoint = CharacteristicPoint.CalculateGasReflect(point, isRight, PbPc, workingGas, innerStreamline)
        if newPoint is None or newPoint.r < 0: #missed the axis
            newPoint = CharacteristicPoint.CaclulateAxisReflect(point, isRight, workingGas)
        innerStreamline = np.append(innerStreamline, newPoint) # append to streamline
        return newPoint, innerStreamline, outerStreamline
    newPoint = CharacteristicPoint.CalculateGasReflect(point, isRight, PambPc, workingGas, outerStreamline)
    outerStreamline = np.append(outerStreamline, newPoint)
    return newPoint, innerStreamline, outerStreamline




def PlotCharacteristicLines(fig: plt.Figure, field: np.ndarray) -> plt.Figure:
    field = np.transpose(field)
    x = np.array([[p.x if p is not None else 0 for p in row] for row in field])
    r = np.array([[p.r if p is not None else 0 for p in row] for row in field])

    ax = fig.axes[0]

    ax.plot(x, r, '-k', linewidth=.5) # L
    ax.grid('on', linestyle='--')

    return fig

def PlotFieldData(fig: plt.Figure, fieldGrid: np.ndarray[CharacteristicPoint], lines: int = 2, stations: int = 5):
    x = np.array([[p.x if p is not None else 0 for p in row] for row in fieldGrid])
    r = np.array([[p.r if p is not None else 0 for p in row] for row in fieldGrid])
    mach = np.array([[p.mach if p is not None else np.nan for p in row] for row in fieldGrid])
    theta = np.array([[p.theta if p is not None else np.nan for p in row] for row in fieldGrid])

    ax = fig.axes[0]

    qx = x[::fieldGrid.shape[0]//stations, ::fieldGrid.shape[1]//lines]
    qy = r[::fieldGrid.shape[0]//stations, ::fieldGrid.shape[1]//lines]

    thetaVx = np.cos(theta[::fieldGrid.shape[0]//stations, ::fieldGrid.shape[1]//lines])
    thetaVy = np.sin(theta[::fieldGrid.shape[0]//stations, ::fieldGrid.shape[1]//lines])

    machContours = ax.contourf(x, r, mach, levels=100, cmap='jet')
    fig.colorbar(machContours, orientation='vertical')

    ax.quiver(qx, qy, thetaVx, thetaVy, scale=25, scale_units='xy', angles='xy', headwidth=3, headlength=5, width=.002, color='black')


def GridifyComplexField(rlines: np.ndarray, llines: np.ndarray) -> np.ndarray:
    R0 = rlines.shape[0]
    L0 = llines.shape[0]

    reflections = (rlines.shape[1] - 1) / (R0 + L0)

    size = int(100 + (R0 + L0)*(reflections + 1)//2)
    X0 = R0
    gridField = np.empty((size, size), dtype=CharacteristicPoint)
    for r, row in enumerate(rlines):
        pos = -1 + X0 - r
        i, j = (0, r + 1)
        rline = True
        for c, point in enumerate(row):
            if gridField[i, j] is None:
                gridField[i, j] = point.clone()
            elif gridField[i, j].terminate or point.terminate:
                gridField[i, j].setTerminate()
            if pos == (R0 + L0):
                rline = not rline
                pos = 0
            i += 1 if rline else 0
            j += 1 if not rline else 0
            pos += 1
        X0 = R0

    X0 = L0
    for r, row in enumerate(llines):
        pos = -1 + X0 - r
        i, j = (r + 1, 0)
        rline = False
        for c, point in enumerate(row):
            gridField[i, j] = point.clone()
            if pos == (R0 + L0):
                rline = not rline
                pos = 0
            i += 1 if rline else 0
            j += 1 if not rline else 0
            pos += 1
    return gridField

def CalculateThrust(exhaust: Gas, Pamb, Tt: Q_, Rt: Q_, Re: Q_, innerStreamline, baseRadius): # TODO trail termaination when can no longer continue
    phi = np.pi/2 + Tt
    Astar = np.pi/np.sin(phi) * (Re**2 - Rt**2)
    
    exitV = gas.MachToVelocity(1, exhaust)
    exitVx = exitV * np.cos(-Tt)
    momThrust = DESIGN.totalmdot * exitVx
    ic(exitV.to(unitReg.feet/unitReg.second))
    ic(momThrust.to(unitReg.pound_force))

    pressThrust = (gas.StagPressRatio(1, exhaust) * exhaust.stagPress - Pamb) * Astar * np.cos(-Tt)
    ic(pressThrust.to(unitReg.pound_force))

    xt: Q_ = (Re - Rt)*np.tan(Tt)

    lasti = 0
    for i in range(len(innerStreamline)):
        if innerStreamline[i].r < baseRadius:
            break
    
    ic(i, len(innerStreamline))
    contPoints = innerStreamline[0:i]

    
    for point in contPoints:
        if type(point.x) == Q_:
            point.x = point.x.magnitude
        if type(point.r) == Q_:
            point.r = point.r.magnitude

    for point in contPoints:
        plt.plot(point.x, point.r, 'xr')
    
    thrusts = []

    for i in range(len(contPoints[:-1])):
        angle = np.pi - np.arctan2(contPoints[i].r - contPoints[i+1].r, contPoints[i].x - contPoints[i+1].x)
        area = Q_(np.pi / np.sin(angle) * (contPoints[i].r**2 - contPoints[i+1].r**2), unitReg.inch**2)
        pressure = gas.StagPressRatio(contPoints[i].mach, exhaust)*exhaust.stagPress
        thrusts.append((pressure - Pamb) * area * np.cos(angle))
    pressureIntegral = sum(thrusts)
    ic(pressureIntegral.to(unitReg.pound_force))

    total = momThrust + pressThrust + pressureIntegral
    ic(total.to(unitReg.pound_force))
    return total