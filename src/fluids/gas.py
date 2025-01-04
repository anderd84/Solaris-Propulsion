import numpy as np
from dataclasses import dataclass
from scipy.optimize import fsolve
from icecream import ic

from general.units import Q_, unitReg
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
    Rgas: Q_
    stagTemp: Q_ = Q_(0, unitReg.degR)
    stagPress: Q_ = Q_(0, unitReg.psi)

    def __init__(self, gamma: float, Rgas: float, **kwargs):
        self.gammaTyp = SpHeatRatio(gamma)
        self.Rgas = Rgas

        if 'T0' in kwargs:
            self.stagTemp = kwargs['T0']
        if 'P0' in kwargs:
            self.stagPress = kwargs['P0']
        if 'Tt' in kwargs:
            self.stagTemp = kwargs['Tt']
        if 'Pt' in kwargs:
            self.stagPress = kwargs['Pt']
    
    def __expr__(self):
        return f"Gas - k:{self.gammaTyp}, R:{self.Rgas}, P0:{self.stagPress}, T0:{self.stagTemp}"

    def getVariableGamma(self, mach) -> SpHeatRatio:
        T = self.stagTemp * (1 + self.gammaTyp[2] * mach**2)
        gammaNext = self.SimpleHarmonicGamma(T)
        gammaPrev = self.gammaTyp.g
        while (abs(gammaNext - gammaPrev) > 1e-6):
            gammaPrev = gammaNext.g
            T = self.stagTemp * (1 + gammaNext[2] * mach**2)
            gammaNext = self.SimpleHarmonicGamma(T)
        return gammaNext
    
    def getChokedArea(self, mdot):
        gamma1 = self.getVariableGamma(1)
        a = mdot*np.sqrt(self.stagTemp)/self.stagPress
        b = np.sqrt(gamma1/self.Rgas)*gamma1[1]**(-gamma1[1]*gamma1[3])
        return a/b

    def SimpleHarmonicGamma(self, temp: float):
        THETA = Q_(5500, unitReg.degR)
        tr = THETA/temp
        return SpHeatRatio(1 + (self.gammaTyp - 1)/(1 + (self.gammaTyp - 1)*((tr)**2)*(np.exp(tr))/(np.exp(tr) - 1)**2))

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

def StagPressRatio(mach, gas: Gas):
    gamma = gas.getVariableGamma(mach)
    return (1 + (gamma[2]*mach**2))**(-gamma[5])

def StagTempRatio(mach, gas: Gas):
    gamma = gas.getVariableGamma(mach)
    return (1 + (gamma[2]*mach**2))**(-1)

def MachToVelocity(mach, gas: Gas):
    gamma = gas.getVariableGamma(mach)
    temp = gas.stagTemp * StagTempRatio(mach, gas)
    ic(temp)
    return mach * np.sqrt(gamma * temp * gas.Rgas)