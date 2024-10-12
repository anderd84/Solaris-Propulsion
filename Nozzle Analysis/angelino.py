import numpy as np
import gas

def CalculateNozzleContour(expansionRatio, radiusBaseNonDim, machExit, gamma):
    # DERIVED VALUES
    angleFlowExit = gas.PrandtlMeyerFunction(machExit, gamma)

    # CALCULATIONS
    mach = np.linspace(1, machExit, 100)
    alpha = angleFlowExit - gas.PrandtlMeyerFunction(mach, gamma) + gas.MachAngle(mach)
    area = gas.Isentropic1DExpansion(1, mach, gamma)

    zeta = (1 - (1 - (area*(1-radiusBaseNonDim**2)*mach*np.sin(alpha)/expansionRatio))**.5) / np.sin(alpha)
    return alpha, zeta