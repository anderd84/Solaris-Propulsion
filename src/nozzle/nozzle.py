import numpy as np
from dataclasses import dataclass

from fluids.gas import MachAngle, SpHeatRatio, PrandtlMeyerFunction, Isentropic1DExpansion

@dataclass
class ContourPoint:
    x: float
    r: float

    def __init__(self, x: float, r: float):
        self.x = x
        self.r = r


def InternalPreExpansion(machE: float, thetaE: float, machLip: float, thetaLip: float, gamma: SpHeatRatio, startPoint: tuple[float, float], expRatio: float, steps: int = 100) -> tuple[np.ndarray, np.ndarray]:
    startX, startR = startPoint
    RrRe = .15
    PhiT = np.pi/2

    Mei = machLip

    print(f"Mei: {Mei}")
    nuei = PrandtlMeyerFunction(Mei, gamma)
    thetaei = thetaLip - nuei
    phiei = thetaei + MachAngle(Mei)

    Rp = np.sqrt(1- (Mei * Isentropic1DExpansion(Mei, gamma) * np.sin(phiei)/expRatio))
    Xp = (Rp - 1)/np.tan(phiei)

    machs = np.linspace(1, Mei, steps)
    nus = PrandtlMeyerFunction(machs, gamma)
    betas = PhiT - np.pi/2 - nus + abs(thetaei)
    L = 2*RrRe*np.sin(.5*betas)

    psi = np.pi - PhiT + nus - (np.pi - betas)/2
    phis = 2*nuei - PrandtlMeyerFunction(machE, gamma) - nus + MachAngle(machs)

    innerR = Rp + L*np.sin(psi)
    innerX = Xp - L*np.cos(psi)

    outerR = np.sqrt(innerR**2 + np.sin(phis)*Isentropic1DExpansion(machs, gamma)*machs/expRatio)
    outerX = innerX + (outerR - innerR)/np.tan(phis)

    return np.array([innerX, innerR]), np.array([outerX, outerR])
    # return Xp, Rp


def RaoContourFormat(contour, scale = 1):
    return np.array([ContourPoint(p.x * scale, p.r * scale) for p in contour[::-1]])