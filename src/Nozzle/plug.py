import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from icecream import ic

from fluids.gas import Gas
import fluids.gas as gas
import Nozzle.rao as rao
from Nozzle import nozzle
from Nozzle import plots
from General.units import Q_, unitReg
import General.design as DESIGN


def CalcPlugLength(machLip: float, theta:float, exhaustGas: Gas, PbPc: float):
    _, _, length = rao.CalculatePlugMetrics(machLip, theta, rao.CalculateMachD(machLip, theta, exhaustGas.gammaTyp, PbPc), exhaustGas.gammaTyp)
    return length

def CreateRaoContour(exhaustGas: Gas, chamberPressure: Q_, designAmbient: Q_, basePress: Q_, lipRadius: Q_, maxSpikeLength: Q_, resolution:int = 50):
    PbPc = basePress / chamberPressure
    ic(PbPc)
    PambPc = designAmbient / chamberPressure

    machLip = fsolve(lambda m: gas.StagPressRatio(m, exhaustGas) - PambPc, 2)[0]
    
    GUESS_T = -.05
    length = CalcPlugLength(machLip, np.deg2rad(GUESS_T), exhaustGas, PbPc)

    if length*lipRadius > maxSpikeLength:
        thetaLip = fsolve(lambda t: CalcPlugLength(machLip, np.deg2rad(t), exhaustGas, PbPc) - maxSpikeLength/lipRadius, -5)[0]
    else:
        thetaLip = GUESS_T

    thetaLip = np.deg2rad(thetaLip)

    areaRatio, Cf, lengthRatio = rao.CalculatePlugMetrics(machLip, thetaLip, rao.CalculateMachD(machLip, thetaLip, exhaustGas.gammaTyp, PbPc), exhaustGas.gammaTyp)

    thetaThroat = rao.CalculateThroatAngle(machLip, thetaLip, 1, exhaustGas.gammaTyp)

    controlSurface: np.ndarray[rao.CharacteristicPoint] = rao.GetControlSurfaceProperties(machLip, thetaLip, lengthRatio, exhaustGas.gammaTyp, resolution)
    expansionFan: np.ndarray[rao.CharacteristicPoint] = rao.GenerateExpansionFan(machLip, 1, thetaThroat, exhaustGas.gammaTyp, resolution)

    field = rao.GenerateFlowField(expansionFan, controlSurface, exhaustGas.gammaTyp)

    radiusThroat = np.sqrt(1 - (1/areaRatio*np.cos(thetaThroat)))

    cont = rao.CalculateContour(field, radiusThroat, thetaThroat)

    field = rao.PruneField(field)

    formatContour = nozzle.RaoContourFormat(cont, lipRadius.to(unitReg.inch).magnitude)

    outputData = {"radiusThroat": Q_(radiusThroat*lipRadius, unitReg.inch), "thetaThroat": Q_(thetaThroat, unitReg.radian), "machLip": machLip, "thetaLip": Q_(thetaLip, unitReg.radian), "areaRatio": areaRatio, "Cf": Cf, "lengthRatio": lengthRatio, "rawContour": cont}

    return formatContour, field, outputData

def GenerateDimPlug(contour: np.ndarray[nozzle.ContourPoint], throatRadius: Q_, throatTheta: Q_, Re: Q_, chamberLength: Q_, baseRadius: Q_, circRes: int = 50):
    designTable = DESIGN.plugDesignTable
    designTable["throatArcRad"] = designTable["throatArcRadFactor"]*Re.to(unitReg.inch).magnitude
    designTable["turnArcRad"] = designTable["turnArcRadFactor"]*Re.to(unitReg.inch).magnitude

    throatRadius = throatRadius.to(unitReg.inch).magnitude
    throatTheta = throatTheta.to(unitReg.radian).magnitude
    Re = Re.to(unitReg.inch).magnitude
    chamberLength = chamberLength.to(unitReg.inch).magnitude
    baseRadius = baseRadius.to(unitReg.inch).magnitude

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

    # convergeLine = np.array([nozzle.ContourPoint(x1CL, r1CL), nozzle.ContourPoint(x2CL, r2CL)])

    #converge arc
    xcCA = x2CL - designTable["turnArcRad"]*np.sin(np.deg2rad(designTable["convergeAngle"]))
    ycCA = r2CL + designTable["turnArcRad"]*np.cos(np.deg2rad(designTable["convergeAngle"]))
    arcArr = np.linspace(np.deg2rad(designTable["convergeAngle"]), 0, circRes)
    convergeArc = np.array([nozzle.ContourPoint(xcCA + designTable["turnArcRad"]*np.sin(a), ycCA - designTable["turnArcRad"]*np.cos(a)) for a in arcArr])

    #straight section
    spiketipx = contour[-1].x
    spiketipr = contour[-1].r
    xSS = np.array([x2CL - chamberLength, x2CL - chamberLength, spiketipx])
    rSS = np.array([convergeArc[-1].r, 0, 0])

    straightSection = np.array([nozzle.ContourPoint(xSS[i], rSS[i]) for i in range(len(xSS))])
    
    fullPlugContour = np.concatenate([contour[::-1], throatArc[1:], convergeArc, straightSection], axis=0)
    fullPlugContour = np.insert(fullPlugContour, 0, nozzle.ContourPoint(spiketipx, 0), axis=0)

    return fullPlugContour

def GenerateDimCowl(throatRadius: Q_, throatTheta: Q_, Re: Q_, chamberLength: Q_, chamberOuter: Q_, thickness: Q_, overchoke: Q_, circRes: int = 50):
    designTable = DESIGN.plugDesignTable
    designTable["throatArcRad"] = designTable["throatArcRadFactor"]*Re.to(unitReg.inch).magnitude
    designTable["turnArcRad"] = designTable["turnArcRadFactor"]*Re.to(unitReg.inch).magnitude

    throatRadius = throatRadius.to(unitReg.inch).magnitude
    throatTheta = throatTheta.to(unitReg.radian).magnitude
    Re = Re.to(unitReg.inch).magnitude
    chamberLength = chamberLength.to(unitReg.inch).magnitude
    chamberOuter = chamberOuter.to(unitReg.inch).magnitude
    thickness = thickness.to(unitReg.inch).magnitude

    absThetaT = abs(throatTheta)
    xt = (Re - throatRadius)*np.tan(throatTheta)
    xc = xt - (np.cos(np.pi/2 - absThetaT)*designTable["throatArcRad"])

    #! New stuff
    thetaL = np.deg2rad(designTable["lipAngle"])
    thetal = np.deg2rad(90 - designTable["straightAngle"]) - absThetaT

    sL = np.sin(thetaL)
    cL = np.cos(thetaL)
    sl = np.sin(thetal)
    cl = np.cos(thetal)
    cotl = 1/np.tan(thetal)
    cscl = 1/np.sin(thetal)

    o = overchoke.to(unitReg.inch).magnitude
    Rinner = chamberOuter
    ic(Rinner)
    Rmax = chamberOuter + thickness
    ic(Rmax)
    ic(xc)

    xm = .05 #! make an input

    Amat = np.array([[0, 0, 1, 0, sL, 0, 0, 0, 0], 
                    [0, 0, 0, 1, -cL, 0, 0, 0 ,0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0], 
                    [0, 0, 0, 0, 0, 1, 1, 0, 0],
                    [0, -cl, 1, 0, -cl, 0, 0, -sl, 0],
                    [1, sl, 0, -1, sl, 0, 0, -cl, 0],
                    [0, 0, 0, 0, 0, 0, -sL, 0, cL],
                    [0, 0, 0, 0, 0, 1, -cL, 0, -sL],
                    [0, 0, cotl, 1, -cscl, 0, 0, 0, 0]])

    bmat = np.array([0, Re, Rinner, Rmax, xc, 0, xm, Re, Re-o*cscl])

    yc, rc, xf, yf, rf, ym, rm, l, L = np.linalg.solve(Amat, bmat)
    #! New stuff

    x1 = 0
    r1 = Re

    arcArr = np.linspace(-thetaL, np.pi/2 - thetal, circRes)
    lipArc = np.array([nozzle.ContourPoint(xf - rf*np.sin(a), yf - rf*np.cos(a)) for a in arcArr])

    plt.plot(0, Re, 'rx')
    
    # lipStraight = np.array([nozzle.ContourPoint(xf - rf*cl, yf - rf*sl), nozzle.ContourPoint(xc + rc*cl, yc + rc*sl)])

    arcArr = np.linspace(thetal, np.pi/2, circRes)
    convergeArc = np.array([nozzle.ContourPoint(xc + rc*np.cos(a), yc + rc*np.sin(a)) for a in arcArr])

    chamber = np.array([nozzle.ContourPoint(xc, yc+rc), nozzle.ContourPoint(xc - chamberLength, yc+rc), nozzle.ContourPoint(xc - chamberLength, Rmax), nozzle.ContourPoint(xm, Rmax)])

    arcArr = np.linspace(np.pi/2, -np.pi/2 + thetaL, circRes)
    manifold = np.array([nozzle.ContourPoint(xm + rm*np.cos(a), ym + rm*np.sin(a)) for a in arcArr])

    final = np.array([nozzle.ContourPoint(0, Re)])


    return np.concatenate([lipArc, convergeArc, chamber, manifold, final], axis=0)