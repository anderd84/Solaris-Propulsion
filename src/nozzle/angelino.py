import numpy as np
import matplotlib.pyplot as plt
from icecream import ic

from fluids import gas
import general.design as DESIGN

def CalculateNozzleContour(expansionRatio, radiusBaseNonDim, machExit, gamma):
    # DERIVED VALUES
    angleFlowExit = gas.PrandtlMeyerFunction(machExit, gamma)

    # CALCULATIONS
    mach = np.linspace(1, machExit, 100)
    alpha = angleFlowExit - gas.PrandtlMeyerFunction(mach, gamma) + gas.MachAngle(mach)
    area = gas.Isentropic1DExpansion(mach, gamma)

    zeta = (1 - (1 - (area*(1-radiusBaseNonDim**2)*mach*np.sin(alpha)/expansionRatio))**.5) / np.sin(alpha)
    return alpha, zeta

def CalculateInternalExternalNozzle(machInlet: float, machExit: float, radiusBaseNonDim: float, exhaust: gas.Gas, steps: int = 100):
    #lip at (0,1)
    gamma = exhaust.gammaTyp
    deltaNu = gas.PrandtlMeyerFunction(machExit, gamma) - gas.PrandtlMeyerFunction(machInlet, gamma)
    thetaInlet = -deltaNu
    thetas = np.linspace(thetaInlet, 0, steps)

    machs = np.linspace(machInlet, machExit, steps)
    alphas = gas.MachAngle(machs)
    angles = thetas - alphas

    expansionRatio = gas.Isentropic1DExpansion(machExit, gamma)
    areas = gas.Isentropic1DExpansion(machs, gamma)

    r = np.sqrt(1 - (machs*areas*np.sin(-angles)/expansionRatio)) # full length, not truncated
    angles = angles[r >= radiusBaseNonDim]
    r = r[r >= radiusBaseNonDim]
    x = (r - 1) / np.tan(angles)



    deltaNu = gas.PrandtlMeyerFunction(machInlet, gamma)
    thetaThroat = thetaInlet + deltaNu
    curveRad = .25
    phiThroat = thetaThroat + np.pi/2
    # thetaInlet =  thetaThroat - deltaNu 
    phiInlet = thetaInlet + gas.MachAngle(machInlet)
    rP = np.sqrt(1 - (machInlet*gas.Isentropic1DExpansion(machInlet, gamma)*np.sin(phiInlet)/expansionRatio))
    xP = (rP - 1) / np.tan(phiInlet)

    machs = np.linspace(1, machInlet, steps)
    nus = gas.PrandtlMeyerFunction(machs, gamma)
    alphas = gas.MachAngle(machs)
    areas = gas.Isentropic1DExpansion(machs, gamma)
    B = thetaThroat - nus + abs(thetaInlet)
    L = 2*curveRad*np.sin(.5*B)
    psi = np.pi - phiThroat + nus - (np.pi - B)/2
    rPlug = rP + L*np.sin(psi)
    xPlug = xP - L*np.cos(psi)

    phis = 2*deltaNu - gas.PrandtlMeyerFunction(machExit, gamma) - nus + alphas
    rCowl = np.sqrt(rPlug**2 + (machs*areas*np.sin(phis)/expansionRatio))
    xCowl = xPlug + (rCowl - rPlug) / np.tan(phis)

    phiThroat = phis[0]
    astarnd = np.pi/np.sin(phiThroat)*(rCowl[0]**2 - rPlug[0]**2)
    scale = np.sqrt(DESIGN.chokeArea.m/astarnd)

    xPlug *= scale
    rPlug *= scale
    xCowl *= scale
    rCowl *= scale
    x *= scale
    r *= scale

    ic(scale)
    ic(np.pi/np.sin(phiThroat)*(rCowl[0]**2 - rPlug[0]**2))
    ic(DESIGN.chokeArea.m)

    for i in range(0, len(x), 10):
        plt.plot([0, x[i]], [scale, r[i]], '-k', linewidth=0.5)
    plt.plot(x, r, '-k')

    plt.plot(np.append(xPlug, x[0]), np.append(rPlug, r[0]), '-r')
    plt.plot(xCowl, rCowl, '-b')

    for i in range(0, len(xCowl), 10):
        plt.plot([xPlug[i], xCowl[i]], [rPlug[i], rCowl[i]], '-k', linewidth=0.5)


    plt.grid(True)
    plt.axis('equal')
    plt.show()