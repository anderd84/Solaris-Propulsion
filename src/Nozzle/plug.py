import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from icecream import ic

from fluids.gas import Gas
import fluids.gas as gas
import Nozzle.rao as rao
from Nozzle import nozzle
from Nozzle import plots


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
    # formatContour = np.append(formatContour, nozzle.ContourPoint((lipRadius - radiusThroat*lipRadius)*np.tan(thetaThroat), 0))
    # formatContour = np.append(formatContour, nozzle.ContourPoint(formatContour[0].x, 0))
    # formatContour = np.append(formatContour, nozzle.ContourPoint(formatContour[0].x, formatContour[0].r))

    outputData = {"radiusThroat": radiusThroat*lipRadius, "thetaThroat": thetaThroat, "machLip": machLip, "thetaLip": thetaLip, "areaRatio": areaRatio, "Cf": Cf, "lengthRatio": lengthRatio, "rawContour": cont}

    return formatContour, field, outputData

def GenerateDimPlug(contour: np.ndarray[nozzle.ContourPoint], throatRadius: float, throatTheta: float, Re: float, chamberLength: float, baseRadius: float, circRes: int = 50):
    designTable = {"throatArcRad": .133*Re, "convergeAngle": 25, "turnArcRad": 2*Re, "straightAngle": 10}

    xt = (Re - throatRadius)*np.tan(throatTheta)
    absThetaT = abs(throatTheta)

    # thoat arc
    ycTA = throatRadius - (np.sin(np.pi/2 - absThetaT)*designTable["throatArcRad"])
    xcTA = xt - (np.cos(np.pi/2 - absThetaT)*designTable["throatArcRad"])
    arcArr = np.linspace(np.pi/2 - absThetaT, np.pi/2 + np.deg2rad(designTable["convergeAngle"]), circRes)
    throatArc = np.array([nozzle.ContourPoint(xcTA + designTable["throatArcRad"]*np.cos(a), ycTA + designTable["throatArcRad"]*np.sin(a)) for a in arcArr])

    # converge line
    x1CL = throatArc[-1].x
    r1CL = throatArc[-1].r

    r2CL = baseRadius + designTable["turnArcRad"]*(1 - np.cos(np.deg2rad(designTable["convergeAngle"])))
    x2CL = x1CL - (r1CL - r2CL)/np.tan(np.deg2rad(designTable["convergeAngle"]))

    convergeLine = np.array([nozzle.ContourPoint(x1CL, r1CL), nozzle.ContourPoint(x2CL, r2CL)])

    #converge arc
    xcCA = x2CL - designTable["turnArcRad"]*np.sin(np.deg2rad(designTable["convergeAngle"]))
    ycCA = r2CL + designTable["turnArcRad"]*np.cos(np.deg2rad(designTable["convergeAngle"]))
    arcArr = np.linspace(np.deg2rad(designTable["convergeAngle"]), 0, circRes)
    convergeArc = np.array([nozzle.ContourPoint(xcCA + designTable["turnArcRad"]*np.sin(a), ycCA - designTable["turnArcRad"]*np.cos(a)) for a in arcArr])

    #straight section
    spiketipx = contour[-1].x
    spiketipr = contour[-1].r
    xSS = np.array([convergeArc[-1].x, xcTA - chamberLength, xcTA - chamberLength, spiketipx, spiketipx])
    rSS = np.array([convergeArc[-1].r, convergeArc[-1].r, 0, 0, spiketipr])

    straightSection = np.array([nozzle.ContourPoint(xSS[i], rSS[i]) for i in range(len(xSS))])
    
    fullPlugContour = np.concatenate([throatArc, convergeLine, convergeArc, straightSection, contour[::-1]], axis=0)

    return fullPlugContour

def GenerateDimCowl(contour: np.ndarray[nozzle.ContourPoint], throatRadius: float, throatTheta: float, Re: float, chamberLength: float, chamberOuter: float, thickness: float, circRes: int = 50):
    designTable = {"throatArcRad": .133*Re, "convergeAngle": 25, "turnArcRad": 2*Re, "straightAngle": 10, "lipAngle": 5}

    absThetaT = abs(throatTheta)
    xt = (Re - throatRadius)*np.tan(throatTheta)
    xc = xt - (np.cos(np.pi/2 - absThetaT)*designTable["throatArcRad"])

    alpha = np.deg2rad(designTable["straightAngle"]) + absThetaT

    Amat = np.array([[0, np.sin(alpha), np.cos(alpha)], 
                     [1, np.cos(alpha), -np.sin(alpha)],
                     [1, 1, 0]])
    Bmat = np.array([-xc, Re, chamberOuter])

    ic(np.linalg.det(Amat))
    yc, r, l = np.linalg.solve(Amat, Bmat)

    x1 = 0
    r1 = Re

    arcArr = np.linspace(alpha, 0, circRes)
    innerArc = np.array([nozzle.ContourPoint(xc + r*np.sin(a), yc + r*np.cos(a)) for a in arcArr])
    xL = innerArc[-1].x - chamberLength
    rL = chamberOuter

    innerCowl = np.insert(innerArc, 0, nozzle.ContourPoint(x1, r1))
    innerCowl = np.append(innerCowl, nozzle.ContourPoint(xL, rL))

    x1 = thickness*np.cos(np.deg2rad(designTable["lipAngle"]))
    r1 = Re + thickness*np.sin(np.deg2rad(designTable["lipAngle"]))

    r = r + thickness
    outerArc = np.array([nozzle.ContourPoint(xc + r*np.sin(a), yc + r*np.cos(a)) for a in arcArr])

    xL = outerArc[-1].x - chamberLength
    rL = chamberOuter + thickness

    outerCowl = np.insert(outerArc, 0, nozzle.ContourPoint(x1, r1))
    outerCowl = np.append(outerCowl, nozzle.ContourPoint(xL, rL))

    fullCowl = np.concatenate([innerCowl, outerCowl[::-1]], axis=0)
    fullCowl = np.append(fullCowl, fullCowl[0])

    return fullCowl