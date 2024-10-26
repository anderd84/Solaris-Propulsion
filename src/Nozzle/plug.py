import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from icecream import ic

from fluids.gas import Gas
import fluids.gas as gas
import Nozzle.rao as rao
from Nozzle import nozzle


def CalcPlugLength(machLip: float, theta:float, exhaustGas: Gas, PbPc: float):
    _, _, length = rao.CalculatePlugMetrics(machLip, theta, rao.CalculateMachD(machLip, theta, exhaustGas.gammaTyp, PbPc), exhaustGas.gammaTyp)
    return length

def CreateRaoContour(exhaustGas: Gas, chamberPressure: float, chamberTemp: float, designAmbient: float, basePress: float, lipRadius: float, maxSpikeLength: float, resolution:int = 50):
    exhaustGas.stagPress = chamberPressure
    exhaustGas.stagTemp = chamberTemp

    PbPc = basePress / chamberPressure
    ic(PbPc)
    PambPc = designAmbient / chamberPressure

    machLip = fsolve(lambda m: gas.StagPressRatio(m, exhaustGas) - PambPc, 2)[0]
    
    GUESS_T = -.05
    length = CalcPlugLength(machLip, np.deg2rad(GUESS_T), exhaustGas, PbPc)

    ic(length*lipRadius)

    if length*lipRadius > maxSpikeLength:
        thetaLip = fsolve(lambda t: CalcPlugLength(machLip, np.deg2rad(t), exhaustGas, PbPc) - maxSpikeLength/lipRadius, -5)[0]
    else:
        thetaLip = GUESS_T

    ic(thetaLip)
    thetaLip = np.deg2rad(thetaLip)

    areaRatio, Cf, lengthRatio = rao.CalculatePlugMetrics(machLip, thetaLip, rao.CalculateMachD(machLip, thetaLip, exhaustGas.gammaTyp, PbPc), exhaustGas.gammaTyp)

    thetaThroat = rao.CalculateThroatAngle(machLip, thetaLip, 1, exhaustGas.gammaTyp)

    controlSurface: np.ndarray[rao.CharacteristicPoint] = rao.GetControlSurfaceProperties(machLip, thetaLip, lengthRatio, exhaustGas.gammaTyp, resolution)
    expansionFan: np.ndarray[rao.CharacteristicPoint] = rao.GenerateExpansionFan(machLip, 1, thetaThroat, exhaustGas.gammaTyp, resolution)

    field = rao.GenerateFlowField(expansionFan, controlSurface, exhaustGas.gammaTyp)

    radiusThroat = np.sqrt(1 - (1/areaRatio*np.cos(thetaThroat)))

    cont = rao.CalculateContour(field, radiusThroat, thetaThroat)

    field = rao.PruneField(field)

    formatContour = nozzle.RaoContourFormat(cont, lipRadius)

    outputData = {"radiusThroat": radiusThroat*lipRadius, "thetaThroat": thetaThroat, "machLip": machLip, "thetaLip": thetaLip, "areaRatio": areaRatio, "Cf": Cf, "lengthRatio": lengthRatio}

    return formatContour, field, outputData


    




