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

    fig = plt.figure()
    
    plt.plot(data[:, :, 2], data[:, :, 0])
    plt.plot(np.transpose(data[:, :, 2]), np.transpose(data[:, :, 0]))
    plt.grid(True)
    plt.show()


def InputDataGenerate(machE, thetaE, gamma, deltaMach, PbPc):
    thetaE = np.deg2rad(thetaE)

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
        LI.eta = 1/(np.tan(LI.alpha)*LI.machStar)
        LI.beta = np.sin(LI.theta)*np.sin(LI.alpha)/(LI.r*np.cos(LI.theta + LI.alpha))
        return LI

    def RightInvarient(self):
        RI = self.clone()
        RI.lambda_ = np.tan(RI.theta - RI.alpha)
        RI.eta = 1/(np.tan(RI.alpha)*RI.machStar)
        RI.beta = np.sin(RI.theta)*np.sin(RI.alpha)/(RI.r*np.cos(RI.theta - RI.alpha))
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
        return CharacteristicPoint(self.x, self.r, self.theta, self.alpha, self.mach, self.machStar, self.lambda_, self.eta, self.beta)
    
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

# tt =  te -ne -nt
# t = tt + m - nt

def CalculateThroatAngle(machE, thetaE, machT, gamma):
    return thetaE - PrandtlMeyerFunction(machE, gamma) + PrandtlMeyerFunction(machT, gamma)

def GenerateExpansionFan(machE, thetaE, machT, thetaT, gamma, arraySize = 100) -> np.ndarray[CharacteristicPoint]:
    machArr = np.linspace(machT, machE, arraySize)
    nuT = PrandtlMeyerFunction(machT, gamma)

    expansionFanArray = np.zeros(arraySize, dtype=CharacteristicPoint)
    for i, mach in enumerate(machArr):
        delta = PrandtlMeyerFunction(mach, gamma) - nuT
        theta = thetaT + delta
        expansionFanArray[i] = CharacteristicPoint(0, 1, theta, MachAngle(mach), mach=mach)

    return expansionFanArray[-2::-1]

def GenerateFlowField(expansionFanArray: np.ndarray[CharacteristicPoint], controlSurfaceArray: np.ndarray[CharacteristicPoint], gamma):
    def CalculateFieldPoint(L: CharacteristicPoint, R: CharacteristicPoint, gamma: float):
        L, R = PrepBoundaryPoints(L, R, gamma)
        N = CharacteristicPoint(0, 0, 0, 0).LRCombine(L, R)

        N.mach = machStar2mach(N.machStar, gamma)
        N.alpha = np.arcsin(1/N.mach)
        for i in range(30):
            NL = N.LeftInvarient()
            NR = N.RightInvarient()

            NLbar = L.AverageInvarients(NL)
            NRbar = R.AverageInvarients(NR)

            NN = CharacteristicPoint(0,0,0,0).LRCombine(NLbar, NRbar)
            NN.mach = machStar2mach(NN.machStar, gamma)
            NN.alpha = np.arcsin(1/NN.mach)

            print(f"percent diff: {np.abs(NN.theta - N.theta)/abs(NN.theta)}")
            if np.abs((NN.theta - N.theta)/(NN.theta)) < 1e-2:
                print(f"converged in {i + 1} iterations")
                return NN
            else:
                N = NN.clone()

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












def CalculateNozzleContour():
    pass