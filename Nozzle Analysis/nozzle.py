import numpy as np
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

    mplot = [machB]
    tplot = [thetaB]

    def isAtD(mach, theta): 
        # test for end
        temp = 2*np.sqrt(mach*mach - 1)/(gamma*mach*mach)
        temp2 = (1+gam2*mach*mach)**(gam5)
        return (np.sin(-2 * theta) - temp + PbPc*temp*temp2) > 0 # want to be positive

    while not isAtD(machA, thetaA):
        # current state calculations
        machA = machB - deltaMach
        tanAlphaA = TanAlpha(machA)
        alphaA = np.arctan(tanAlphaA)
        machStarA = MachStar(machA)

        mach2 = machStarA**2
        gammaTemp = gam4*mach2/(mach2 - 1)
        temp1 = (1-gam6*mach2)/(mach2-1)
        temp2 = (-2*A)/machStarA*np.sqrt(temp1)
        temp3 = (4*A*A)/mach2*temp1 - 4*gammaTemp*((A*A/mach2)-1)
        temp3 = np.sqrt(temp3)
        temp4 = 2*gammaTemp

        sinTheta = (temp2+temp3)/temp4
        thetaA = np.arcsin(sinTheta)

        rRatioA = B / ((machStarA*np.sin(-thetaA))**2 * tanAlphaA*(1+gam2*machA**2)**(-gam3))

        t1 = (1/((1+gam2*machB**2)*gam4)**gam3)
        t2 = np.sqrt((gam1*machB*machB)/(1+gam2*machB*machB))*rRatioA
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

    return areaRatio, Cf, lengthRatio, machA, np.rad2deg(thetaA), mplot, tplot



