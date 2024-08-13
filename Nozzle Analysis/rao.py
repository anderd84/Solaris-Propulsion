from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from gas import mach2machStar, machStar2mach, PrandtlMeyerFunction, MachAngle

def GenerateInputMatrix(machArray, thetaArray, gamma, deltaMach, PbPc):
    data = np.zeros((len(machArray), len(thetaArray), 3))
    for i, mach in enumerate(machArray):
        for j, theta in enumerate(thetaArray):
            print(mach, theta)
            data[i, j, :] = InputDataGenerate(mach, theta, gamma, deltaMach, PbPc)
    return data

def GenerateInputChart(machArray, thetaArray, gamma, deltaMach, PbPc):
    data = GenerateInputMatrix(machArray, thetaArray, gamma, deltaMach, PbPc)

    plt.plot(data[:, :, 2], data[:, :, 0])
    plt.plot(np.transpose(data[:, :, 2]), np.transpose(data[:, :, 0]))
    plt.grid(True)
    plt.show()


def InputDataGenerate(machE, thetaE, gamma, deltaMach, PbPc):
    # useful gamma terms
    gam1 = (gamma + 1) / 2
    gam2 = (gamma - 1) / 2
    gam3 = 1 / (gamma - 1)
    gam4 = 2 / (gamma + 1)
    gam5 = gamma / (gamma - 1)
    gam6 = (gamma - 1) / (gamma + 1)

    TanAlpha = lambda mach: 1 / np.sqrt(mach**2 - 1)
    cot = lambda x: np.cos(x) / np.sin(x)

    machStarE = mach2machStar(machE, gamma) 
    tanAlphaE = TanAlpha(machE)
    alphaE = np.arctan(tanAlphaE)

    A = machStarE * (np.cos(thetaE) + tanAlphaE*np.sin(-thetaE))
    B = (machStarE * np.sin(-thetaE))**2 * (tanAlphaE*(1+gam2*machE**2)**(-gam3))

    # need to track 4 states while integrating, mach, theta, alpha, and rRatio
    # B terms represent the state from the previous iteration, A terms are for the current iteration
    machA = 1
    thetaA = 0
    alphaA = 0
    rRatioA = 1

    machB = machE
    thetaB = thetaE
    alphaB = alphaE
    rRatioB = 1

    sumAreaRatio = 0
    sumLengthRatio = 0
    sumForce = 0

    areaRatio = 0
    Cf = 0
    lengthRatio = 0

    mplot = []
    tplot = []
    aplot = []
    rplot = []
    lplot = []

    def isAtD(mach, theta):
        # test for end
        test = (2/(gamma*mach*mach))*np.sqrt(mach*mach - 1) - np.sin(-2*theta)
        c = (1 + gam2*mach*mach)**gam5
        c1 = (2/(gamma*mach*mach))*np.sqrt(mach*mach - 1)
        c2 = PbPc * c1 * c
        return (test - c2 < 0) # want to be negative

    while not isAtD(machA, thetaA):
        # current state calculations
        machA = machB - deltaMach
        tanAlphaA = TanAlpha(machA)
        alphaA = np.arctan(tanAlphaA)
        machStarA = mach2machStar(machA, gamma)

        mach2 = machStarA**2
        gammaTemp = gam4*mach2/(mach2 - 1)
        temp1 = (1-gam6*mach2)/(mach2 - 1)
        temp2 = (-2*A)/machStarA*np.sqrt(temp1)
        temp3 = np.sqrt((4*A*A)/mach2*temp1 - 4*gammaTemp*((A*A/mach2)-1))
        temp4 = 2*gammaTemp

        sinThetaA = (temp2+temp3)/temp4
        thetaA = np.arcsin(sinThetaA)

        rRatioA = B / ((machStarA * np.sin(-thetaA))**2 * (tanAlphaA*(1+gam2*machA**2)**(-gam3)))

        t1 = (1/((1+gam2*machB**2)*gam4)**gam3)
        t2 = np.sqrt((gam1*machB*machB)/(1+gam2*machB*machB))*rRatioB
        t2 = t2*np.sin(alphaB)/np.sin(thetaB-alphaB)

        t3 = (1/((1+gam2*machA**2)*gam4)**gam3)
        t4 = np.sqrt((gam1*machA*machA)/(1+gam2*machA*machA))*rRatioA
        t4 = t4*np.sin(alphaA)/np.sin(thetaA-alphaA)

        dRR = rRatioA - rRatioB

        sumAreaRatio += (t1*t2 + t3*t4)*dRR
        sumLengthRatio += .5*(cot(thetaB - alphaB) + cot(thetaA - alphaA))*dRR

        f1 = (1+gam2*machB**2)**(-gam5)
        f2 = (1 + gamma*machB**2*((np.sin(alphaB)*np.cos(-thetaB))/(np.sin(-thetaB+alphaB))))*rRatioB

        f3 = (1+gam2*machA**2)**(-gam5)
        f4 = (1 + gamma*machA**2*((np.sin(alphaA)*np.cos(-thetaA))/(np.sin(-thetaA+alphaA))))*rRatioA

        sumForce += (f1*f2 + f3*f4)*(-dRR)

        areaRatio = 1/sumAreaRatio
        Cf = areaRatio*sumForce
        lengthRatio = sumLengthRatio

        machB = machA
        thetaB = thetaA
        alphaB = alphaA
        rRatioB = rRatioA
        
        mplot.append(machB)
        tplot.append(thetaB)
        aplot.append(areaRatio)
        rplot.append(rRatioB)
        lplot.append(lengthRatio)

    return areaRatio, Cf, lengthRatio, aplot, rplot, lplot, tplot, mplot

@dataclass
class CharacteristicPoint:
    x: float
    r: float
    theta: float
    alpha: float
    mach: float = 0
    machStar: float = 0
    lambda_: float = 0
    eta: float = 0
    beta: float = 0

    def LeftInvarient(self):
        LI = self.clone()
        LI.lambda_ = np.tan(LI.theta + LI.alpha)
        LI.eta = 1/np.tan(LI.alpha)/LI.machStar
        LI.beta = np.sin(LI.theta)*np.sin(LI.alpha)/(LI.r*np.cos(LI.theta + LI.alpha))
        return LI

    def RightInvarient(self):
        RI = self.clone()
        RI.lambda_ = np.tan(RI.theta - RI.alpha)
        RI.eta = 1/np.tan(RI.alpha)/RI.machStar
        RI.beta = np.sin(RI.theta)*np.sin(RI.alpha)/(RI.r*np.sin(RI.theta - RI.alpha))
        return RI

    def AverageInvarients(self, N: 'CharacteristicPoint'):
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

    def clone(self):
        return CharacteristicPoint(x=self.x, r=self.r, theta=self.theta, alpha=self.alpha, mach=self.mach, machStar=self.machStar, lambda_=self.lambda_, eta=self.eta, beta=self.beta)
    
def GetControlSurfaceProperties(machE, thetaE, lengthRatio, gamma, arraySize = 100):
    gam1 = (gamma + 1) / 2
    gam2 = (gamma - 1) / 2
    gam3 = 1 / (gamma - 1)
    gam4 = 2 / (gamma + 1)
    gam5 = gamma / (gamma - 1)
    gam6 = (gamma - 1) / (gamma + 1)

    TanAlpha = lambda mach: 1 / np.sqrt(mach**2 - 1)

    machStarE = mach2machStar(machE, gamma)
    tanAlphaE = TanAlpha(machE)

    A = machStarE * (np.cos(thetaE) + tanAlphaE*np.sin(-thetaE))
    B = (machStarE * np.sin(-thetaE))**2 * (tanAlphaE*(1+gam2*machE**2)**(-gam3))

    # start at x=0, go down to x=lenghRatio
    xArr = np.linspace(0, lengthRatio, arraySize)

    controlSurfaceArray = np.zeros(arraySize, dtype=CharacteristicPoint)
    controlSurfaceArray[0] = CharacteristicPoint(0, 1, thetaE, np.arctan(tanAlphaE), mach=machE)

    def CalcTheta(M):
        machStarA = mach2machStar(M, gamma)

        mach2 = machStarA**2
        gammaTemp = gam4*mach2/(mach2 - 1)
        temp1 = (1-gam6*mach2)/(mach2 - 1)
        temp2 = (-2*A)/machStarA*np.sqrt(temp1)
        temp3 = np.sqrt(np.abs((4*A*A)/mach2*temp1 - 4*gammaTemp*((A*A/mach2)-1)))
        temp4 = 2*gammaTemp

        sinThetaA = (temp2+temp3)/temp4
        theta = np.arcsin(sinThetaA)
        return theta
    
    def CalcB(mach, theta, rRatio):
        machStar = mach2machStar(mach, gamma)
        tanAlpha = TanAlpha(mach)
        return (machStar * np.sin(-theta))**2 * (tanAlpha*(1+gam2*mach**2)**(-gam3)) * rRatio

    for i, x in enumerate(xArr[1:]):
        dx = x - xArr[i]
        # dr = np.tan(thetaArr[i] - alphaArr[i]) * dx
        dr = np.tan(controlSurfaceArray[i].theta - controlSurfaceArray[i].alpha) * dx
        # r = rArr[i] + dr
        r = controlSurfaceArray[i].r + dr
        mach = fsolve(lambda M: CalcB(M, CalcTheta(M), r) - B, controlSurfaceArray[i].mach, xtol=1e-12)[0]
        theta = CalcTheta(mach)
        alpha = np.arcsin(1/mach)

        controlSurfaceArray[i+1] = CharacteristicPoint(x, r, theta, alpha, mach=mach)
        
    return controlSurfaceArray

def CalculateThroatAngle(machE, thetaE, machT, gamma):
    return thetaE - PrandtlMeyerFunction(machE, gamma)# + PrandtlMeyerFunction(machT, gamma)

def GenerateExpansionFan(machE, thetaE, machT, thetaT, gamma, arraySize = 100) -> np.ndarray[CharacteristicPoint]:
    machArr = np.linspace(machT, machE, arraySize)
    nuT = PrandtlMeyerFunction(machT, gamma)

    expansionFanArray = np.zeros(arraySize, dtype=CharacteristicPoint)
    for i, mach in enumerate(machArr):
        delta = PrandtlMeyerFunction(mach, gamma)# - nuT
        theta = thetaT + delta
        expansionFanArray[i] = CharacteristicPoint(0, 1, theta, MachAngle(mach), mach=mach)

    return expansionFanArray[-2::-1]

def GenerateFlowField(expansionFanArray: np.ndarray[CharacteristicPoint], controlSurfaceArray: np.ndarray[CharacteristicPoint], gamma):
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

            print(f"percent diff: {np.abs(NN.theta - N.theta)/abs(NN.theta)}")
            if np.abs((NN.theta - N.theta)/(NN.theta)) < 1e-4:
                print(f"converged in {i + 1} iterations")
                NN.mach = machStar2mach(NN.machStar, gamma)
                NN.alpha = MachAngle(NN.mach)
                return NN
            else:
                N = NN.clone()
                del NN, NL, NR, NLbar, NRbar

        print("did not converge")
        return N

    def PrepBoundaryPoints(L: CharacteristicPoint, R: CharacteristicPoint, gamma: float):
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
            print(f"calculating at point {i}, {j}")
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

def CalculateContour(field: np.ndarray[CharacteristicPoint], Rt: float, Tt: float):
    pt1 = field[0, -2]
    pt2 = field[1, -1]
    L = field[0, -1]

    cont = [L]

    m = 0
    n = 0

    rows = field.shape[0]

    def eqn(pt1: CharacteristicPoint, pt2: CharacteristicPoint, L: CharacteristicPoint, ax: CharacteristicPoint):
        return (L.r - pt1.r + (pt1.x - ax)*(pt1.r-pt2.r)/(pt1.x - pt2.x))/(L.x - ax) - np.tan((pt1.x-ax)/(pt1.x-pt2.x)*pt2.theta+(ax-pt2.x)/(pt1.x-pt2.x)*pt1.theta)

    def calcNext(pt1, pt2, L):
        ax2 = fsolve(lambda ax: eqn(pt1, pt2, L, ax), pt2.x)[0]
        ar2 = pt1.r - (pt1.x - ax2)*(pt1.r - pt2.r)/(pt1.x - pt2.x)
        print(ax2, ar2)
        return CharacteristicPoint(ax2, ar2, 0, 0)

    while m < rows - 2:
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

    return cont

def PruneUnderContour(field, contour):
    xtol = field[0,-1].x/field.shape[0]
    rtol = field[0,0].r/field.shape[1]
    print(xtol, rtol)
    for i in range(1, field.shape[0]):
        for j in range(1, field.shape[1]):
            bad = False
            for p in contour:
                if abs(p.x - field[i, j].x) < xtol and p.r - field[i, j].r < rtol:
                    bad = True
            if not bad:
                field[i, j] = field[i, j-1].clone()
    return field

def distance(p1, p2):
    return np.sqrt((p1.x - p2.x)**2 + (p1.r - p2.r)**2)