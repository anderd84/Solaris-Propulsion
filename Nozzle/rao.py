from dataclasses import dataclass
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.integrate as integrate
from gas import mach2machStar, machStar2mach, PrandtlMeyerFunction, MachAngle, SpHeatRatio

def cot(x: float) -> float:
    return np.cos(x) / np.sin(x)

def CalcTheta(mach: float, gamma: SpHeatRatio, A: float) -> float:
    machStar = mach2machStar(mach, gamma)
    mach2 = machStar**2
    gammaTemp = gamma[4]*mach2/(mach2 - 1)
    temp1 = (1-gamma[6]*mach2)/(mach2 - 1)
    temp2 = (-2*A)/machStar*np.sqrt(temp1)
    temp3 = np.sqrt(abs((4*A*A)/mach2*temp1 - 4*gammaTemp*((A*A/mach2)-1)))
    temp4 = 2*gammaTemp

    sinTheta = (temp2+temp3)/temp4
    return np.arcsin(sinTheta)

def GenerateInputMatrix(machArray: npt.ArrayLike, thetaArray: npt.ArrayLike, gamma: SpHeatRatio, PbPc: float) -> np.ndarray:
    data = np.zeros((len(machArray), len(thetaArray), 3))
    for i, mach in enumerate(machArray):
        for j, theta in enumerate(thetaArray):
            data[i, j, :] = CalculatePlugMetrics(mach, theta, CalculateMachD(mach, theta, gamma, PbPc), gamma, steps=10000)
    return data

def GenerateInputChart(data: np.ndarray) -> None:
    plt.contourf(np.transpose(data[:, :, 2]), np.transpose(data[:, :, 0]), data[:, :, 1]/np.max(data[:,:,1]), levels=10, cmap='jet')
    plt.plot(data[:, :, 2], data[:, :, 0], 'k', linewidth=.5)
    plt.plot(np.transpose(data[:, :, 2]), np.transpose(data[:, :, 0]), 'k', linewidth=.5)
    plt.colorbar()
    plt.grid(True)
    plt.show()

def CalculateMachD(machE: float, thetaE: float, gamma: SpHeatRatio, PbPc: float, tol: float = 1e-8) -> float:
    machStarE = mach2machStar(machE, gamma) 
    alphaE = MachAngle(machE)
    tanAlphaE = np.tan(alphaE)

    A = machStarE * (np.cos(thetaE) + tanAlphaE*np.sin(-thetaE))
    B = (machStarE * np.sin(-thetaE))**2 * (tanAlphaE*(1+gamma[2]*machE**2)**(-gamma[3]))

    def EndCondition(mach):
        theta = CalcTheta(mach, gamma, A)
        test = (2/(gamma*mach*mach))*np.sqrt(abs(mach*mach - 1)) - np.sin(-2*theta)
        c = (1 + gamma[2]*mach*mach)**gamma[5]
        c1 = (2/(gamma*mach*mach))*np.sqrt(abs(mach*mach - 1))
        c2 = PbPc * c1 * c
        return test - c2
    
    machD = fsolve(EndCondition, machE - .2, xtol=tol)[0]

    return machD

def CalculatePlugMetrics(machE: float, thetaE: float, machD: float, gamma: SpHeatRatio, steps: int = 100) -> tuple:    
    machStarE = mach2machStar(machE, gamma) 
    alphaE = MachAngle(machE)
    tanAlphaE = np.tan(alphaE)

    A = machStarE * (np.cos(thetaE) + tanAlphaE*np.sin(-thetaE))
    B = (machStarE * np.sin(-thetaE))**2 * (tanAlphaE*(1+gamma[2]*machE**2)**(-gamma[3]))

    machs = np.linspace(machD, machE, steps)
    machStars = mach2machStar(machs, gamma)
    thetas = CalcTheta(machs, gamma, A)
    alphas = MachAngle(machs)

    rRatios = B / ((machStars * np.sin(-thetas))**2 * (np.tan(MachAngle(machs))*(1+gamma[2]*machs**2)**(-gamma[3])))

    t1 = -(1/((1+gamma[2]*machs**2)*gamma[4])**gamma[3])
    t2 = np.sqrt((gamma[1]*machs*machs)/(1+gamma[2]*machs*machs))*rRatios
    t2 *= np.sin(alphas)/np.sin(thetas-alphas)

    f1 = (1+gamma[2]*machs**2)**(-gamma[5])
    f2 = (1 + gamma*machs**2*((np.sin(alphas)*np.cos(-thetas))/(np.sin(-thetas+alphas))))*rRatios

    areaRatio = 1/integrate.trapezoid(2*t1*t2, rRatios)
    lengthRatio = integrate.trapezoid(cot(-thetas + alphas), rRatios)
    Cf = areaRatio*integrate.trapezoid(2*f1*f2, rRatios)

    return areaRatio, Cf, lengthRatio

@dataclass
class CharacteristicPoint:
    """
    Stores data relating to a point on a characteristic line
    Also has some functions for calculating the left and right invarients <- idk if thats what they're called

    ### Attributes:
    1. x: float
        x value of the point
    2. r: float
        r value of the point
    3. theta: float
        flow angle at the point
    4. alpha: float
        mach angle at the point
    5. mach: float
        mach number at the point
    6. machStar: float
        mach star at the point

    7. lambda_: float
        invarient variable
    8. eta: float  
        invarient variable
    9. beta: float
        invarient variable

    ### Methods:
    1. LeftInvarient
        calculates the left invarient variables
    2. RightInvarient
        calculates the right invarient variables
    3. AverageInvarients
        calculates the average invarient variables between two points for a new point
    4. LRCombine
        use invarient variables from a left and right point to approximate states at the point between them
    5. clone
        returns a copy of the object
    """
    x: float
    r: float
    theta: float
    alpha: float
    mach: float = 0
    machStar: float = 0
    lambda_: float = 0
    eta: float = 0
    beta: float = 0

    def LeftInvarient(self) -> 'CharacteristicPoint':
        LI = self.clone()
        LI.lambda_ = np.tan(LI.theta + LI.alpha)
        LI.eta = 1/np.tan(LI.alpha)/LI.machStar
        LI.beta = np.sin(LI.theta)*np.sin(LI.alpha)/(LI.r*np.cos(LI.theta + LI.alpha))
        return LI

    def RightInvarient(self) -> 'CharacteristicPoint':
        RI = self.clone()
        RI.lambda_ = np.tan(RI.theta - RI.alpha)
        RI.eta = 1/np.tan(RI.alpha)/RI.machStar
        RI.beta = np.sin(RI.theta)*np.sin(RI.alpha)/(RI.r*np.sin(RI.theta - RI.alpha))
        return RI

    def AverageInvarients(self, N: 'CharacteristicPoint') -> 'CharacteristicPoint':
        lambda_ = (self.lambda_ + N.lambda_) / 2
        eta = (self.eta + N.eta) / 2
        beta = (self.beta + N.beta) / 2

        return CharacteristicPoint(x=self.x, r=self.r, theta=self.theta, alpha=self.alpha, mach=self.mach, machStar=self.machStar, lambda_=lambda_, eta=eta, beta=beta)
    
    def LRCombine(self, L: 'CharacteristicPoint', R: 'CharacteristicPoint') -> 'CharacteristicPoint':
        self.x = ((R.lambda_*R.x - L.lambda_*L.x) + (L.r - R.r)) / (R.lambda_ - L.lambda_)
        self.r = L.r - L.lambda_*(L.x - self.x)
        self.machStar = (R.theta - L.theta + L.eta*L.machStar + R.eta*R.machStar - R.beta*(R.r - self.r) - L.beta*(L.x - self.x)) / (L.eta + R.eta)
        self.theta = L.theta - L.eta*(L.machStar - self.machStar) + L.beta*(L.x - self.x)
        return self

    def clone(self) -> 'CharacteristicPoint':
        return CharacteristicPoint(x=self.x, r=self.r, theta=self.theta, alpha=self.alpha, mach=self.mach, machStar=self.machStar, lambda_=self.lambda_, eta=self.eta, beta=self.beta)
    
def GetControlSurfaceProperties(machE: float, thetaE: float, lengthRatio: float, gamma: SpHeatRatio, arraySize: int = 100) -> np.ndarray[CharacteristicPoint]:
    """
    this function also needs to be finalized, equations are same so error propogates
    """
    TanAlpha = lambda mach: 1 / np.sqrt(mach**2 - 1)

    machStarE = mach2machStar(machE, gamma)
    tanAlphaE = TanAlpha(machE)

    A = machStarE * (np.cos(thetaE) + tanAlphaE*np.sin(-thetaE))
    B = (machStarE * np.sin(-thetaE))**2 * (tanAlphaE*(1+gamma[2]*machE**2)**(-gamma[3]))

    # start at x=0, go down to x=lenghRatio
    xArr = np.linspace(0, lengthRatio, arraySize)

    controlSurfaceArray = np.zeros(arraySize, dtype=CharacteristicPoint)
    controlSurfaceArray[0] = CharacteristicPoint(0, 1, thetaE, np.arctan(tanAlphaE), mach=machE)
    
    def CalcB(mach: float, theta: float, rRatio: float) -> float:
        machStar = mach2machStar(mach, gamma)
        tanAlpha = TanAlpha(mach)
        return (machStar * np.sin(-theta))**2 * (tanAlpha*(1+gamma[2]*mach**2)**(-gamma[3])) * rRatio

    for i, x in enumerate(xArr[1:]):
        dx = x - xArr[i]
        # dr = np.tan(thetaArr[i] - alphaArr[i]) * dx
        dr = np.tan(controlSurfaceArray[i].theta - controlSurfaceArray[i].alpha) * dx
        # r = rArr[i] + dr
        r = controlSurfaceArray[i].r + dr
        mach = fsolve(lambda M: CalcB(M, CalcTheta(M, gamma, A), r) - B, controlSurfaceArray[i].mach, xtol=1e-12)[0]
        theta = CalcTheta(mach, gamma, A)
        alpha = np.arcsin(1/mach)

        controlSurfaceArray[i+1] = CharacteristicPoint(x, r, theta, alpha, mach=mach)
        
    return controlSurfaceArray

def CalculateThroatAngle(machE: float, thetaE: float, machT: float, gamma: SpHeatRatio) -> float:
    return thetaE - PrandtlMeyerFunction(machE, gamma) + PrandtlMeyerFunction(machT, gamma)

def GenerateExpansionFan(machE: float, machT: float, thetaT: float, gamma: SpHeatRatio, arraySize: int = 100) -> np.ndarray[CharacteristicPoint]: # look into making the linspace angles instead of mach
    machArr = np.linspace(machT, machE, arraySize)
    nuT = PrandtlMeyerFunction(machT, gamma)

    expansionFanArray = np.zeros(arraySize, dtype=CharacteristicPoint)
    for i, mach in enumerate(machArr):
        delta = PrandtlMeyerFunction(mach, gamma) - nuT
        theta = thetaT + delta
        expansionFanArray[i] = CharacteristicPoint(0, 1, theta, MachAngle(mach), mach=mach)

    return expansionFanArray[-2::-1]

def GenerateFlowField(expansionFanArray: np.ndarray[CharacteristicPoint], controlSurfaceArray: np.ndarray[CharacteristicPoint], gamma: SpHeatRatio) -> np.ndarray[CharacteristicPoint]:
    def CalculateFieldPoint(L: CharacteristicPoint, R: CharacteristicPoint, gamma: float):
        L, R = PrepBoundaryPoints(L, R, gamma)
        N = CharacteristicPoint(0, 0, 0, 0).LRCombine(L, R)

        for i in range(30):
            N.mach = machStar2mach(N.machStar, gamma)
            N.alpha = MachAngle(N.mach)
            
            NL = N.LeftInvarient()
            NR = N.RightInvarient()

            NLbar = L.AverageInvarients(NL)
            NRbar = R.AverageInvarients(NR)

            NN = CharacteristicPoint(0,0,0,0).LRCombine(NLbar, NRbar)

            # print(f"percent diff: {np.abs(NN.theta - N.theta)/abs(NN.theta)}")
            if np.abs((NN.theta - N.theta)/(NN.theta)) < 1e-4:
                # print(f"converged in {i + 1} iterations")
                NN.mach = machStar2mach(NN.machStar, gamma)
                NN.alpha = MachAngle(NN.mach)
                return NN
            else:
                N = NN.clone()
                del NN, NL, NR, NLbar, NRbar

        # print("did not converge")
        return N

    def PrepBoundaryPoints(L: CharacteristicPoint, R: CharacteristicPoint, gamma: SpHeatRatio) -> tuple[CharacteristicPoint, CharacteristicPoint]:
        L.machStar = mach2machStar(L.mach, gamma)
        L = L.LeftInvarient()
        R.machStar = mach2machStar(R.mach, gamma)
        R = R.RightInvarient()
        return L, R

    field = np.empty((len(expansionFanArray) + 1, len(controlSurfaceArray)), dtype=CharacteristicPoint) # v axis is expansion fan, h axis is control surface
    field[0, :] = controlSurfaceArray
    field[1:, 0] = expansionFanArray
    
    for j in range(1, len(controlSurfaceArray)):
        for i in range(1, len(expansionFanArray) + 1):
            # print(f"calculating at point {i}, {j}")
            field[i, j] = CalculateFieldPoint(field[i-1, j], field[i, j-1], gamma)
            # field[i, j] = CharacteristicPoint(0,0,0,0)
    return field

def PruneField(field: np.ndarray[CharacteristicPoint]) -> np.ndarray[CharacteristicPoint]:
    for i in range(1, field.shape[0]):
        for j in range(1, field.shape[1]):
            if field[i, j].r > field[i, j-1].r:
                field[i, j] = field[i, j-1].clone()
            if np.isnan(field[i, j].r):
                field[i, j] = field[i, j-1].clone()
    return field

def CalculateContour(field: np.ndarray[CharacteristicPoint], Rt: float, Tt: float) -> np.ndarray[CharacteristicPoint]:
    pt1 = field[0, -2]
    pt2 = field[1, -1]
    L = field[0, -1]

    cont = [L]

    m = 0
    n = 0

    rows, cols = field.shape


    def eqn(pt1: CharacteristicPoint, pt2: CharacteristicPoint, L: CharacteristicPoint, ax: CharacteristicPoint) -> float:
        return (L.r - pt1.r + (pt1.x - ax)*(pt1.r-pt2.r)/(pt1.x - pt2.x))/(L.x - ax) - np.tan((pt1.x-ax)/(pt1.x-pt2.x)*pt2.theta+(ax-pt2.x)/(pt1.x-pt2.x)*pt1.theta)

    def calcNext(pt1: CharacteristicPoint, pt2: CharacteristicPoint, L: CharacteristicPoint) -> CharacteristicPoint:
        ax2 = fsolve(lambda ax: eqn(pt1, pt2, L, ax), pt2.x)[0]
        ar2 = pt1.r - (pt1.x - ax2)*(pt1.r - pt2.r)/(pt1.x - pt2.x)
        # print(ax2, ar2)
        return CharacteristicPoint(ax2, ar2, np.arctan((ar2 - L.r)/(ax2 - L.x)), 0)

    while m < rows - 2 and n > -(cols - 2):
        pt1 = field[0 + m, -2 + n]
        pt2 = field[1 + m, -1 + n]
        try:
            Lnew = calcNext(pt1, pt2, cont[-1]) # 1
        except:
            print(f"Failed at {(m, n)}")
            break
        if Lnew.x > max(pt2.x, pt1.x) or Lnew.x < min(pt2.x, pt1.x):
            m += 1
            n += 1
        else:
            cont.append(Lnew)
            n -= 1

    cont.append(CharacteristicPoint((1 - Rt)*np.tan(Tt), Rt, 0, 0))

    return np.array(cont)

def PruneUnderContour(field: np.ndarray[CharacteristicPoint], contour: np.ndarray[CharacteristicPoint]) -> np.ndarray[CharacteristicPoint]:
    xtol = field[0,-1].x/field.shape[0]
    rtol = field[0,0].r/field.shape[1]
    # print(xtol, rtol)
    for i in range(1, field.shape[0]):
        for j in range(1, field.shape[1]):
            bad = False
            for p in contour:
                if abs(p.x - field[i, j].x) < xtol and p.r - field[i, j].r < rtol:
                    bad = True
            if not bad:
                field[i, j] = field[i, j-1].clone()
    return field

def distance(p1: CharacteristicPoint, p2: CharacteristicPoint) -> float:
    return np.sqrt((p1.x - p2.x)**2 + (p1.r - p2.r)**2)