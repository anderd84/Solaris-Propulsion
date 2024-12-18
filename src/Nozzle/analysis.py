import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.optimize import fsolve
import fluids.gas as gas
from fluids.gas import MachAngle, SpHeatRatio, mach2machStar, machStar2mach, Gas
import Nozzle.nozzle as nozzle
import matrix_viewer as mv
from icecream import ic
import Nozzle.plots as plots
import Nozzle.config as config

def CalculateSimpleField(contour, PambPc, PbPc, gamma, Mt, Tt, steps = 100, reflections = 2):
    phit = np.pi/2 + Tt
    Me = np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2])
    machs = np.linspace(Mt, Me, steps)
    temps = 1/(1 + gamma[2]*machs**2)
    press = temps**(gamma[5])
    charLineAngles = Tt + (gas.PrandtlMeyerFunction(machs, gamma) - gas.PrandtlMeyerFunction(Mt, gamma)) - gas.MachAngle(machs)
    thetaEnd = Tt + (gas.PrandtlMeyerFunction(Me, gamma) - gas.PrandtlMeyerFunction(Mt, gamma))

    nu1 = gas.PrandtlMeyerFunction(Mt, gamma)
    nu2 = gas.PrandtlMeyerFunction(Me, gamma)

    print(f"Me: {Me}, Mt: {Mt}")
    print(f"ThetaEnd: {np.rad2deg(thetaEnd)}, Tt: {np.rad2deg(Tt)}")

    endpoints, bounds = CalculateEndpoints(contour, charLineAngles)
    # gas.weakObliqueShock(2, np.deg2rad(10), gamma)
    temps = temps[bounds[0]:bounds[1]]
    machs = machs[bounds[0]:bounds[1]]
    # dA = CalculateDiffArea(endpoints)


    plt.figure()
    plt.plot([p.x for p in endpoints], [p.r for p in endpoints], '.r')
    plt.plot([p.x for p in contour], [p.r for p in contour], '--k', linewidth=0.5)
    # xgrid, ygrid, chord = SimpleCharacteristicsMesh(endpoints, machs)
    # plt.contourf(xgrid, ygrid, chord, levels=100, cmap='jet')

    ColorLine(endpoints[:], temps[:], 'jet', plt.gca())
    plt.plot([1, 1 + np.cos(thetaEnd)], [1, 1 + np.sin(thetaEnd)], '-r', linewidth=2)
    plt.plot([0, np.sin(Tt)], [1, 1 - np.cos(Tt)], '-b', linewidth=2)

    plt.plot([0, 10*np.cos(Tt + nu1)], [1, 1+10*np.sin(Tt + nu1)], '--r', linewidth=1)
    plt.plot([0, 10*np.cos(Tt + nu2 - nu1)], [1, 1+10*np.sin(Tt + nu2 - nu1)], '--g', linewidth=1)

    for ep in endpoints:
        plt.plot([0, ep.x], [1, ep.r], '--k', linewidth=0.1)
    


    plt.xlim(contour[0].x, contour[-1].x)
    plt.ylim(0, 2)
    plt.gca().set_aspect('equal')
    plt.show()

def CalculateEndpoints(contour, angles):
    rayLength = np.sqrt((0 - contour[-1].x)**2 + (1 - contour[-1].r)**2)
    
    endpoints = np.empty(len(angles), dtype=nozzle.ContourPoint)

    lowerBound = 0
    for i, angle in enumerate(angles):
        for j in range(len(contour) - 1):
            a, b, c, d = contour[j].x, contour[j].r, contour[j+1].x, contour[j+1].r
            rayEndx = 0 + rayLength * np.cos(angle)
            rayEndr = 1 + rayLength * np.sin(angle)
            x = ((d - (b - d)/(a - c)*c) - 1)/((rayEndr - 1)/(rayEndx) - (b - d)/(a - c))
            if a <= x <= c:
                endpoints[i] = nozzle.ContourPoint(x, (rayEndr - 1)/(rayEndx)*x + 1)
                break
            if j == len(contour) - 2: 
                if endpoints[i - 1] is not None:
                    return endpoints[lowerBound:i], (lowerBound, i)
                else:
                    lowerBound = i + 1
    return endpoints, (lowerBound, len(endpoints))

def ColorLine(points, colorFunc, colormap, axes):
    x, y = [p.x for p in points], [p.r for p in points]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=colormap, norm=plt.Normalize(min(colorFunc), max(colorFunc)), array=colorFunc)
    line = axes.add_collection(lc)
    plt.colorbar(line, ax=axes)
    return line

def SimpleCharacteristicsMesh(contour, chordFunc):
    xgrid = np.array([p.x for p in contour])
    xgrid = np.array([np.zeros(len(contour)), xgrid])

    rgrid = np.array([p.r for p in contour])
    rgrid = np.array([np.ones(len(contour)), rgrid])

    chord = np.array([chordFunc[:len(contour)], chordFunc[:len(contour)]])

    return xgrid, rgrid, chord

def CalculateDiffArea(contour):
    dA = np.zeros(len(contour) - 1)
    for i in range(len(contour) - 1):
        a, b, c, d = contour[i].x, contour[i].r, contour[i+1].x, contour[i+1].r
        theta = np.arctan((d - b)/(c - a))
        dA[i] = abs(np.pi/np.sin(theta)*(b**2 - d**2))
    return dA





@dataclass
class CharacteristicPoint:
    x: float
    r: float
    theta: float
    machStar: float
    s: float = 0
    mach: float = 0
    alpha: float = 0

    F: float = 0
    G: float = 0
    H: float = 0
    J: float = 0

    def __repr__(self) -> str:
        return f"CP @ ({self.x:.3f}, {self.r:.3f}), theta: {np.rad2deg(self.theta):.3f}, s: {self.s:.3f}, mach: {self.mach:.3f}"

    def clone(self):
        return CharacteristicPoint(self.x, self.r, self.theta, self.machStar, self.s, self.mach, self.alpha)

    def CalculateRightVariant(self, gas: Gas) -> 'CharacteristicPoint': # I characteristic
        Rgas = gas.Rgas
        gamma = gas.gammaTyp
        RV = self.clone()
        RV.F = np.tan(RV.theta - RV.alpha)
        RV.G = 1/np.tan(RV.alpha)/RV.machStar
        RV.H = -np.sin(RV.theta)*np.sin(RV.alpha)/(RV.r*np.sin(RV.theta - RV.alpha))
        RV.J = np.sin(RV.alpha)*np.cos(RV.alpha)/(Rgas * gamma)
        return RV

    def CalculateLeftVariant(self, gas: Gas) -> 'CharacteristicPoint': # II characteristic
        Rgas = gas.Rgas
        gamma = gas.gammaTyp
        LV = self.clone()
        LV.F = np.tan(LV.theta + LV.alpha)
        LV.G = -1/np.tan(LV.alpha)/LV.machStar
        LV.H = np.sin(LV.theta)*np.sin(LV.alpha)/(LV.r*np.cos(LV.theta + LV.alpha))
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
            if np.abs((NN.theta - N.theta)/(NN.theta)) < 1e-6:
                # ic(f"Converged in {i} iterations")
                return NN
            else:
                N = NN.clone()
                del NN, NL, NR, L2, R2

        # ic("Failed to converge")
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

        Amat = np.array([[1, R.G], [1, L.G]])
        b = np.array([[R.theta + R.G*R.machStar - R.H*(r - R.r) - R.J*(s - R.s)], [L.theta + L.G*L.machStar - L.H*(x - L.x) - L.J*(s - L.s)]])
        X = np.linalg.solve(Amat, b)
        theta = X[0, 0]
        machStar = X[1, 0]
        mach = machStar2mach(machStar, workingGas.gammaTyp)
        alpha = MachAngle(mach)

        return CharacteristicPoint(x, r, theta, machStar, s, mach, alpha)

    @staticmethod #TODO Make this iterative
    def CalculateSolidReflect(point: 'CharacteristicPoint', isRight: bool, contour: np.ndarray[nozzle.ContourPoint], PbPc: float, workingGas: Gas) -> 'CharacteristicPoint':
        ic(isRight)
        PV = point.CalculateRightVariant(workingGas) if isRight else point.CalculateLeftVariant(workingGas)
        # ic(PV)
        intersect, theta = CharacteristicPoint.CalculateSolidBoundaryIntersect(PV, contour, isRight)
        ic(intersect)
        if intersect is None:
            return PV.clone()#CharacteristicPoint.CalculateGasReflect(point, isRight, PbPc, workingGas)
        
        x = intersect[0]
        r = intersect[1]
        s = 0
        machStar = PV.machStar + (-(theta - PV.theta) - PV.H*(r - PV.r) - PV.J*(s - PV.s))/PV.G

        return CharacteristicPoint(x, r, theta, machStar, s, machStar2mach(machStar, workingGas.gammaTyp), MachAngle(machStar2mach(machStar, workingGas.gammaTyp)))

    @staticmethod #TODO make this iterative
    def CalculateGasReflect(point: 'CharacteristicPoint', isRight, PambPc, workingGas: Gas, streamline: np.ndarray['CharacteristicPoint']) -> 'CharacteristicPoint':
        PV = point.CalculateRightVariant(workingGas) if isRight else point.CalculateLeftVariant(workingGas)
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
        theta = PV.theta - PV.G*(machStar - PV.machStar) - PV.H*(x - PV.x) - PV.J*(s - PV.s)

        return CharacteristicPoint(x, r, theta, machStar, s, machStar2mach(machStar, workingGas.gammaTyp), MachAngle(machStar2mach(machStar, workingGas.gammaTyp)))



    @staticmethod
    def CaclulateBaseReflect():
        pass

    @staticmethod
    def CalculateSolidBoundaryIntersect(point: 'CharacteristicPoint', contour: np.ndarray[nozzle.ContourPoint], isRight: bool) -> tuple[float, float] | None:        
        Sx, Sy = point.x, point.r
        angle = point.theta - point.alpha if isRight else point.theta + point.alpha
        ic(np.rad2deg(angle))

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

def CalculateComplexField(contour, PambPc, PbPc, workingGas: Gas, Mt, Tt, scale = 1, Rsteps = 20, Lsteps = 0, reflections = 3):
    gamma = workingGas.gammaTyp
    Me = np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2])
    thetaExit = Tt + gas.PrandtlMeyerFunction(Me, gamma) - gas.PrandtlMeyerFunction(Mt, gamma)

    outerStreamLine = np.array([CharacteristicPoint(0, scale, thetaExit, mach2machStar(Me, gamma), 0, Me, MachAngle(Me))])

    rLines = np.empty((Rsteps, 1 + (Lsteps + Rsteps)*reflections), dtype=CharacteristicPoint)
    lLines = np.empty((Lsteps, 1 + (Lsteps + Rsteps)*reflections), dtype=CharacteristicPoint)

    rLines[:, 0] = np.transpose(GenerateExpansionFan(Me, Mt, Tt, workingGas, Rsteps, scale))
    lLines[:, 0] = np.transpose(GenerateExpansionFan(Me, Mt, Tt, workingGas, Lsteps, scale))

    ic(rLines)

    for i in range(reflections):
        PropogateRegionAll(rLines, lLines, workingGas, i)
        rlines, lines, outerStreamLine = ReflectionRegionAll(rLines, lLines, contour, PambPc, PbPc, outerStreamLine, workingGas, i)

    return np.concatenate((rLines, lLines), axis=0), outerStreamLine

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

def ReflectionRegionAll(rLines: np.ndarray[CharacteristicPoint], lLines: np.ndarray[CharacteristicPoint], contour, PambPc, PbPc, streamline, workingGas: Gas, reflection: int):
    R0: int = rLines.shape[0]
    L0: int = lLines.shape[0]
    
    rlines, streamline = ReflectionRegion(rLines, R0, L0, contour, PambPc, PbPc, streamline, workingGas, reflection, True)
    llines, streamline = ReflectionRegion(lLines, L0, R0, contour, PambPc, PbPc, streamline, workingGas, reflection, False)
    return rlines, llines, streamline

def ReflectionRegion(lines: np.ndarray[CharacteristicPoint], X0, Y0, contour, PambPc, PbPc, streamline, workingGas: Gas, reflection, startAsRight: bool): # startAsRight is true if region is being calculated in rlines
    start = 1 + Y0 + (reflection)*(X0 + Y0) # reflection should start at 0
    region = np.empty((X0+1, X0+1), dtype=CharacteristicPoint)
    region[1:,0] = lines[:,start-1]
    region[0,1:] = np.transpose(lines[:,start-1])

    # plt.ion()
    # fig = post.CreateNonDimPlot()
    # post.PlotContour(fig, contour, 0, 0)

    for i in range(1, X0+1):
        for j in range(1, i+1):
            isRight = not (reflection % 2 == 0) ^ startAsRight
            if i == j:
                reflectOrigin = region[i, j-1] # previous point in the same line
                region[i,j] = CharacteristicPoint.CalculateSolidReflect(reflectOrigin, isRight, contour, PbPc, workingGas) if isRight else CharacteristicPoint.CalculateGasReflect(reflectOrigin, isRight, PambPc, workingGas, streamline)
                if not isRight:
                    streamline = np.append(streamline, region[i,j])
            else:
                region[i,j] = CharacteristicPoint.CalculateFieldPoint(region[i-1,j], region[i,j-1], workingGas) if isRight else CharacteristicPoint.CalculateFieldPoint(region[i,j-1], region[i-1,j], workingGas)
                region[j, i] = region[i, j].clone()

            # if i != j:
            #     fig.axes[0].plot(region[i-1, j].x, region[i-1, j].r, 'or')
            #     fig.axes[0].plot(region[i, j-1].x, region[i, j-1].r, 'ob')
            #     fig.axes[0].plot(region[i, j].x, region[i, j].r, 'og')
            #     fig.axes[0].legend(['', 'Nozzle', 'Left', 'Right', 'Current'])
            # else:
            #     fig.axes[0].plot(region[i, j-1].x, region[i, j-1].r, 'ob')
            #     fig.axes[0].plot(region[i, j].x, region[i, j].r, 'og')
            #     fig.axes[0].legend(['', 'Nozzle', 'Previous', 'Current'])

            # PlotComplexField(fig, region)
            # fig.canvas.draw()
            # fig.canvas.flush_events()
            # plt.waitforbuttonpress()
            # xlims = fig.axes[0].get_xlim()
            # ylims = fig.axes[0].get_ylim()
            # fig.axes[0].clear()
            # post.PlotContour(fig, contour, 0, 0)
            # fig.axes[0].set_xlim(xlims)
            # fig.axes[0].set_ylim(ylims)
    lines[:, start:start+X0] = region[1:,1:]
    return lines, streamline

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

    # ax.quiver(qx, qy, thetaVx, thetaVy, scale=25, scale_units='xy', angles='xy', headwidth=3, headlength=5, width=.002, color='black')


def GridifyComplexField(rlines: np.ndarray, llines: np.ndarray) -> np.ndarray:
    R0 = rlines.shape[0]
    L0 = llines.shape[0]

    reflections = (rlines.shape[1] - 1) / (R0 + L0)

    size = int(1 + (R0 + L0)*(reflections + 1)//2)
    X0 = R0
    gridField = np.empty((size, size), dtype=CharacteristicPoint)
    for r, row in enumerate(rlines):
        pos = -1 + X0 - r
        i, j = (0, r + 1)
        rline = True
        for c, point in enumerate(row):
            gridField[i, j] = point.clone()
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
