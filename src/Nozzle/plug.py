import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from icecream import ic

from fluids.gas import Gas
import fluids.gas as gas
import Nozzle.rao as rao



def CalcPlugLength(machLip: float, theta:float, exhaustGas: Gas, PbPc: float):
    _, _, length = rao.CalculatePlugMetrics(machLip, theta, rao.CalculateMachD(machLip, theta, exhaustGas.gammaTyp, PbPc), exhaustGas.gammaTyp)
    return length

def CreateRaoPlug(exhaustGas: Gas, chamberPressure: float, chamberTemp: float, designAmbient: float, basePress: float, lipRadius: float, maxSpikeLength: float):
    exhaustGas.P0 = chamberPressure
    exhaustGas.T0 = chamberTemp
    machLip = fsolve(lambda m: gas.StagPressRatio(m, exhaustGas), 2)
    PbPc = basePress / chamberPressure
    PambPc = designAmbient / chamberPressure
    
    length = CalcPlugLength(machLip, np.deg2rad(0), exhaustGas)

    if length*lipRadius > maxSpikeLength:
        thetaE = fsolve(lambda t: CalcPlugLength(machLip, np.deg2rad(t), exhaustGas, PbPc)*lipRadius - maxSpikeLength, -15)
        ic(thetaE)
    




