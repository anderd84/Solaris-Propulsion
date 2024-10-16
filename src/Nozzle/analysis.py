import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.optimize import fsolve
import src.fluids.gas as gas
from src.fluids.gas import MachAngle, SpHeatRatio, machStar2mach, Gas
import Nozzle.nozzle as nozzle

def CalculateSimpleField(contour, PambPc, PbPc, gamma, Mt, Tt, steps = 100, reflections = 3):
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

    def clone(self):
        return CharacteristicPoint(self.x, self.r, self.theta, self.machStar, self.s, self.mach, self.alpha)

    def CalculateRightVariant(self, gas: Gas): # I characteristic
        Rgas = gas.Rgas
        gamma = gas.gammaTyp
        RV = self.clone()
        RV.F = np.tan(RV.theta - RV.alpha)
        RV.G = 1/np.tan(RV.alpha)/RV.machStar
        RV.H = -np.sin(RV.theta)*np.sin(RV.alpha)/(RV.r*np.sin(RV.theta - RV.alpha))
        RV.J = np.sin(RV.alpha)*np.cos(RV.alpha)/(Rgas * gamma)
        return RV

    def CalculateLeftVariant(self, gas: Gas): # II characteristic
        Rgas = gas.Rgas
        gamma = gas.gammaTyp
        LV = self.clone()
        LV.F = np.tan(LV.theta + LV.alpha)
        LV.G = -1/np.tan(LV.alpha)/LV.machStar
        LV.H = np.sin(LV.theta)*np.sin(LV.alpha)/(LV.r*np.cos(LV.theta + LV.alpha))
        LV.J = -np.sin(LV.alpha)*np.cos(LV.alpha)/(Rgas * gamma)
        return LV
    
    def setCoefficients(self, F, G, H, J):
        self.F = F
        self.G = G
        self.H = H
        self.J = J

    def NextIterationPoint(self, N: 'CharacteristicPoint'):
        newF = (self.F + N.F)/2
        newG = (self.G + N.G)/2
        newH = (self.H + N.H)/2
        newJ = (self.J + N.J)/2
        nextPoint = self.clone()
        nextPoint.setCoefficients(newF, newG, newH, newJ)

        return nextPoint

    @staticmethod
    def CalculateFieldPoint(L: 'CharacteristicPoint', R: 'CharacteristicPoint', gas: Gas): # TODO overhaul for new coefficient calcs
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
        mach = machStar2mach(machStar, gas.gammaTyp)
        alpha = MachAngle(mach)

        return CharacteristicPoint(x, r, theta, machStar, s, mach, alpha)

    @staticmethod
    def CalculateSolidReflect(point: 'CharacteristicPoint', isRight: bool, contour: np.ndarray[nozzle.ContourPoint], PbPc: float, gas: Gas):
        PV = point.CalculateRightVariant(gas) if isRight else point.CalculateLeftVariant(gas)
        intersect, theta = CharacteristicPoint.CalculateSolidBoundaryIntersect(point, contour)
        if intersect is None:
            return CharacteristicPoint.CalculateGasReflect(point, isRight, PbPc, gas)
        
        x = intersect[0]
        r = intersect[1]
        s = 0
        machStar = PV.machStar + (-(theta - PV.theta) - PV.H*(r - PV.r) - PV.J*(s - PV.s))/PV.G

        return CharacteristicPoint(x, r, theta, machStar, s, machStar2mach(machStar, gas.gammaTyp), MachAngle(machStar2mach(machStar)))

    @staticmethod
    def CalculateGasReflect(point, isRight, PambPc, gas):
        pass

    @staticmethod
    def CalculateSolidBoundaryIntersect(point: 'CharacteristicPoint', contour: np.ndarray[nozzle.ContourPoint]) -> tuple[float, float] | None:
        rayLength = np.sqrt((0 - contour[-1].x)**2 + (1 - contour[-1].r)**2)
        
        Bx, By = 0, 0
        angle = np.atan(point.F)

        for j in range(len(contour) - 1):
            a, b, c, d = contour[j].x, contour[j].r, contour[j+1].x, contour[j+1].r
            rayEndx = 0 + rayLength * np.cos(angle)
            rayEndr = 1 + rayLength * np.sin(angle)
            x = ((d - (b - d)/(a - c)*c) - 1)/((rayEndr - 1)/(rayEndx) - (b - d)/(a - c))
            if a <= x <= c:
                Bx = x
                By = (rayEndr - 1)/(rayEndx)*x + 1
                return (Bx, By), np.tan((c - a)/(d - b))
        return None


def CalculateComplexField(contour, PambPc, PbPc, gamma, Mt, Tt, Rsteps = 10, Lsteps = 0, reflections = 3):
    Me = np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2])

    rLines = np.empty((Rsteps, 1 + (Lsteps + Rsteps)*reflections), dtype=CharacteristicPoint)
    lLines = np.empty((Lsteps, 1 + (Lsteps + Rsteps)*reflections), dtype=CharacteristicPoint)

    rLines[:, 0] = np.transpose(GenerateExpansionFan(Me, Mt, Tt, gamma, Rsteps))
    lLines[:, 0] = np.transpose(GenerateExpansionFan(Me, Mt, Tt, gamma, Lsteps))

    for i in range(reflections):
        PropogateRegionAll(rLines, lLines, gamma, i)
        ReflectionRegionAll(rLines, lLines, contour, PambPc, PbPc, gamma, i)

    return np.concatenate((rLines, lLines), axis=1)

def GenerateExpansionFan(machE: float, machT: float, thetaT: float, gamma: SpHeatRatio, arraySize: int):
    machs = np.linspace(machT, machE, arraySize)

    thetas = thetaT + gas.PrandtlMeyerFunction(machs, gamma) - gas.PrandtlMeyerFunction(machT, gamma)

    expansionFanArray = np.empty(arraySize, dtype=CharacteristicPoint)
    for i, mach in enumerate(machs):
        expansionFanArray[i] = CharacteristicPoint(0, 1, mach, thetas[i])

    return expansionFanArray

def PropogateRegionAll(rLines: np.ndarray[CharacteristicPoint], lLines: np.ndarray[CharacteristicPoint], gamma: gas.SpHeatRatio, reflection: int):
    R0: int = np.size(rLines)[0]
    L0: int = np.size(lLines)[0]

    start = 1 + (reflection)*(R0 + L0) # reflection should start at 0
    region = np.empty((R0+1, L0+1), dtype=CharacteristicPoint)
    region[1:,0] = rLines[:,start-1]
    region[0,1:] = np.transpose(lLines[:,start-1])
    for i in range(1, R0+1):
        for j in range(1, L0+1):
            region[i,j] = CharacteristicPoint.CalculateNewPoint(region[i-1,j], region[i,j-1], gamma) if reflection % 2 == 0 else CharacteristicPoint.CalculateNewPoint(region[i,j-1], region[i-1,j], gamma)
    
    rLines[:,start:start+L0] = region[1:,1:]
    lLines[:,start:start+R0] = np.transpose(region[1:,1:])

def ReflectionRegionAll(rLines: np.ndarray[CharacteristicPoint], lLines: np.ndarray[CharacteristicPoint], contour, PambPc, PbPc, gamma: gas.SpHeatRatio, reflection: int):
    R0: int = np.size(rLines)[0]
    L0: int = np.size(lLines)[0]

    ReflectionRegion(rLines, R0, L0, contour, PambPc, PbPc, gamma, reflection, True)
    ReflectionRegion(lLines, L0, R0, contour, PambPc, PbPc, gamma, reflection, False)

def ReflectionRegion(lines: np.ndarray[CharacteristicPoint], X0, Y0, contour, PambPc, PbPc, gamma, reflection, startAsRight: bool): # startAsRight is true if region is being calculated in rlines
    start = 1 + Y0 + (reflection)*(X0 + Y0) # reflection should start at 0
    region = np.empty((X0+1, X0+1), dtype=CharacteristicPoint)
    region[1:,0] = lines[:,start-1]
    region[0,1:] = np.transpose(lines[:,start-1])
    for i in range(1, X0+1):
        for j in range(1, i+1):
            if i == j:
                region[i,j] = CharacteristicPoint.CalculateSolidReflect(region[i-1,j], region[i,j-1], contour, PbPc, gamma) if ((reflection % 2 == 0) ^ (not startAsRight)) else CharacteristicPoint.CalculateGasReflect(region[i,j-1], region[i-1,j], PambPc, gamma)
            else:
                region[i,j] = CharacteristicPoint.CalculateNewPoint(region[i-1,j], region[i,j-1], gamma) if reflection % 2 == 0 else CharacteristicPoint.CalculateNewPoint(region[i,j-1], region[i-1,j], gamma)
    
    lines[:,start:start+X0] = region[1:,1:]



