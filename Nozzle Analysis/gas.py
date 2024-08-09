import numpy as np

def PrandtlMeyerFunction(M, gamma):
    a = np.sqrt((gamma+1)/(gamma-1))
    b = np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1)))
    c = np.arctan(np.sqrt(M**2-1))
    return a*b-c

def MachAngle(M):
    return np.arcsin(1/M)

def Isentropic1DExpansion(Astar, M, gamma):
    return Astar / M * (1 + (gamma-1)/2*M**2)**((gamma+1)/(2*(gamma-1)))

def mach2machStar(mach, gamma):
    return np.sqrt((((gamma+1)/2)*mach**2) / (1 + ((gamma-1)/2)*mach**2))

def machStar2mach(machStar, gamma):
    return np.sqrt(((2/(gamma+1))*machStar**2)/(1-((gamma-1)/(gamma+1))*machStar**2))