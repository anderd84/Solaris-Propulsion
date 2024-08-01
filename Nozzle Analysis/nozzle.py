import numpy as np
import matplotlib.pyplot as plt
import gas

def Angelino(expansionRatio, radiusBaseNonDim, machExit, gamma):
    # DERIVED VALUES
    angleFlowExit = gas.PrandtlMeyerFunction(machExit, gamma)

    # CALCULATIONS
    mach = np.linspace(1, machExit, 100)
    alpha = angleFlowExit - gas.PrandtlMeyerFunction(mach, gamma) + gas.MachAngle(mach)
    area = gas.Isentropic1DExpansion(1, mach, gamma)

    zeta = (1 - (1 - (area*(1-radiusBaseNonDim**2)*mach*np.sin(alpha)/expansionRatio))**.5) / np.sin(alpha)
    return alpha, zeta

# steps : 
# 1. find max angle and create linspace between inital angle and it
# 2. 
def Rao(expansionRatio, length, gamma):
    pass

def RaoGenerateInputMatrix(machArray, thetaArray, gamma, deltaMach, PbPc):
    data = np.zeros((len(machArray), len(thetaArray), 3))
    for i, mach in enumerate(machArray):
        for j, theta in enumerate(thetaArray):
            print(mach, theta)
            data[i, j, :] = RaoInputDataGenerate(mach, theta, gamma, deltaMach, PbPc)
    return data

def RaoGenerateInputChart(machArray, thetaArray, gamma, deltaMach, PbPc):
    data = RaoGenerateInputMatrix(machArray, thetaArray, gamma, deltaMach, PbPc)

    fig = plt.figure()
    
    plt.plot(data[:, :, 2], data[:, :, 0])
    plt.plot(np.transpose(data[:, :, 2]), np.transpose(data[:, :, 0]))
    plt.grid(True)
    plt.show()


def RaoInputDataGenerate(machE, thetaE, gamma, deltaMach, PbPc):
    thetaE = np.deg2rad(thetaE)

    # useful gamma terms
    gam1 = (gamma + 1) / 2
    gam2 = (gamma - 1) / 2
    gam3 = 1 / (gamma - 1)
    gam4 = 2 / (gamma + 1)
    gam5 = gamma / (gamma - 1)
    gam6 = (gamma - 1) / (gamma + 1)

    MachStar = lambda mach: np.sqrt((gam1*mach**2) / (1 + gam2*mach**2))
    TanAlpha = lambda mach: 1 / np.sqrt(mach**2 - 1)
    cot = lambda x: np.cos(x) / np.sin(x)

    machStarE = MachStar(machE)
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
        c = (1 + gam2 * mach * mach)**gam5
        c1 = (2/(gamma*mach*mach))*np.sqrt(mach*mach - 1)
        c2 = PbPc * c1 * c
        return (test - c2 < 0) # want to be positive

    while not isAtD(machA, thetaA):
        # current state calculations
        machA = machB - deltaMach
        tanAlphaA = TanAlpha(machA)
        alphaA = np.arctan(tanAlphaA)
        machStarA = MachStar(machA)

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



