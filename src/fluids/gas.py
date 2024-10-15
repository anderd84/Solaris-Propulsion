import numpy as np
from dataclasses import dataclass
from scipy.optimize import fsolve

@dataclass
class SpHeatRatio:
    def __init__(self, gamma):
        self.g = gamma
        self.calcShorthands()

    def calcShorthands(self):
        self.s = [0 for i in range(7)]
        self.s[1] = (self.g + 1) / 2
        self.s[2] = (self.g - 1) / 2
        self.s[3] = 1 / (self.g - 1)
        self.s[4] = 2 / (self.g + 1)
        self.s[5] = self.g / (self.g - 1)
        self.s[6] = (self.g - 1) / (self.g + 1)

    def __getitem__(self, key):
        return self.s[key]

    def __add__(self, other):
        return self.g + other
    
    def __sub__(self, other):
        return self.g - other
    
    def __mul__(self, other):
        return self.g * other
    
    def __truediv__(self, other):
        return self.g / other
    
    def __pow__(self, other):
        return self.g ** other
    
    def __rmul__(self, other):
        return other * self.g
    
    def __neg__(self):
        return -self.g

class Gas:
    gammaTyp: SpHeatRatio
    Rgas: float
    T0: float = 0

    def __init__(self, gamma: float, Rgas: float):
        self.gammaTyp = SpHeatRatio(gamma)
        self.Rgas = Rgas

    def getVariableGamma(self, mach):
        return self.gammaTyp

def PrandtlMeyerFunction(M, gamma):
    a = np.sqrt((gamma+1)/(gamma-1))
    b = np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1)))
    c = np.arctan(np.sqrt(M**2-1))
    return a*b-c

def MachAngle(M):
    return np.arcsin(1/M)

def Isentropic1DExpansion(M, gamma: SpHeatRatio):
    return 1/M*((1 + gamma[2]*M**2)/gamma[1])**(gamma[1]*gamma[3])

def mach2machStar(mach, gamma: SpHeatRatio):
    return np.sqrt(abs((((gamma+1)/2)*mach**2) / (1 + ((gamma-1)/2)*mach**2)))

def machStar2mach(machStar, gamma: SpHeatRatio):
    return np.sqrt(abs(((2/(gamma+1))*machStar**2)/(1-((gamma-1)/(gamma+1))*machStar**2)))

def obliqueShock(mach, delta, gamma: SpHeatRatio) -> tuple[float, float, float]:
    betaWeak = fsolve(lambda b: np.tan(delta) - (2/np.tan(b)*(((mach * np.sin(b))**2 - 1)/((mach**2 * (gamma + np.cos(2*b))) + 2))), delta)[0]
    betaStrong = fsolve(lambda b: np.tan(delta) - (2/np.tan(b)*(((mach * np.sin(b))**2 - 1)/((mach**2 * (gamma + np.cos(2*b))) + 2))), np.pi/2)[0]

    mach2 = np.sqrt(((gamma - 1)*(mach*np.sin(betaWeak))**2 + 2)/(2*gamma*(mach*np.sin(betaWeak))**2 - (gamma - 1)))/np.sin(betaWeak - delta)
    return betaWeak, betaStrong, mach2