import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

Me = 4.688
Te = np.deg2rad(-5.4)
gamma = 1.4
PbPc = 0

gam1 = (gamma + 1) / 2
gam2 = (gamma - 1) / 2
gam3 = 1 / (gamma - 1)
gam4 = 2 / (gamma + 1)
gam5 = gamma / (gamma - 1)
gam6 = (gamma - 1) / (gamma + 1)

MachStar = lambda mach: np.sqrt(np.abs((gam1*mach**2) / (1 + gam2*mach**2)))
TanAlpha = lambda mach: 1 / np.sqrt(mach**2 - 1)

machStarE = MachStar(Me)
tanAlphaE = TanAlpha(Me)
alphaE = np.arctan(tanAlphaE)

A = machStarE * (np.cos(Te) - tanAlphaE*np.sin(Te))

def EndFunc(M) -> float:
    theta = CalcTheta(M)

    return DCondition(M, theta)

def DCondition(mach, theta):
    test = (2/(gamma*mach*mach))*np.sqrt(np.abs(mach*mach - 1)) - np.sin(-2*theta)
    c = (1 + gam2 * mach * mach)**gam5
    c1 = (2/(gamma*mach*mach))*np.sqrt(np.abs(mach*mach - 1))
    c2 = PbPc * c1 * c
    return test - c2

def CalcTheta(M):
    machStarA = MachStar(M)

    mach2 = machStarA**2
    gammaTemp = gam4*mach2/(mach2 - 1)
    temp1 = (1-gam6*mach2)/(mach2 - 1)
    temp2 = (-2*A)/machStarA*np.sqrt(temp1)
    temp3 = np.sqrt(np.abs((4*A*A)/mach2*temp1 - 4*gammaTemp*((A*A/mach2)-1)))
    temp4 = 2*gammaTemp

    sinThetaA = (temp2+temp3)/temp4
    theta = np.arcsin(sinThetaA)
    return theta

func = lambda M: CalcTheta(M) - np.arctan(-.154)

val = fsolve(func, 4.5)[0]
print(val)

iv = np.linspace(4, Me+.5, 1000)
dv = np.rad2deg(CalcTheta(iv))

plt.plot(iv, dv)
plt.axvline(Me)
plt.grid(True)
plt.show()

cond = DCondition(val, CalcTheta(val))

print(cond)



# end = fsolve(EndFunc, [4.5], xtol=1e-15)[0]
# print(end)
# print(np.tan(CalcTheta(end)))