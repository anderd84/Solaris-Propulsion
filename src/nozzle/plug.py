import numpy as np
from scipy.optimize import fsolve
from icecream import ic

from fluids.gas import Gas
from fluids import gas
from general.units import Q_, unitReg
import general.design as DESIGN
from nozzle import rao
from nozzle import nozzle

def CalcPlugLength(machLip: float, theta:float, exhaustGas: Gas, PbPc: float):
    _, _, length = rao.CalculatePlugMetrics(machLip, theta, rao.CalculateMachD(machLip, theta, exhaustGas.gammaTyp, PbPc), exhaustGas.gammaTyp)
    return length

def CreateRaoContour(exhaustGas: Gas, chamberPressure: Q_, designAmbient: Q_, basePress: Q_, lipRadiusGuess: Q_, maxSpikeLength: Q_, resolution:int = 50):
    PbPc = basePress / chamberPressure
    PambPc = designAmbient / chamberPressure

    machLip = fsolve(lambda m: gas.StagPressRatio(m, exhaustGas) - PambPc, 2)[0]
    ic(machLip)
    
    GUESS_T = -.1
    length = CalcPlugLength(machLip, np.deg2rad(GUESS_T), exhaustGas, PbPc)

    ic(length)
    ic(maxSpikeLength.magnitude)

    if length > maxSpikeLength.magnitude:
        thetaLip = fsolve(lambda t: CalcPlugLength(machLip, np.deg2rad(t), exhaustGas, PbPc) - maxSpikeLength.magnitude, -3)[0]
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

    phi = np.pi/2 + thetaThroat
    lipRadChoke = np.sqrt(DESIGN.chokeArea * np.sin(phi) / (np.pi * (1 - radiusThroat**2)))

    formatContour = nozzle.RaoContourFormat(cont, lipRadChoke.to(unitReg.inch).magnitude)

    outputData = {"radiusLip": Q_(lipRadChoke, unitReg.inch), "radiusThroat": Q_(radiusThroat*lipRadChoke, unitReg.inch), "thetaThroat": Q_(thetaThroat, unitReg.radian), "machLip": machLip, "thetaLip": Q_(thetaLip, unitReg.radian), "areaRatio": areaRatio, "Cf": Cf, "lengthRatio": lengthRatio, "rawContour": cont}

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
    coolingFloor = DESIGN.coolingChannelWallDist.to(unitReg.inch).magnitude
    coolingHeight = DESIGN.coolingChannelHeightPlug.to(unitReg.inch).magnitude
    thickness = DESIGN.plugThickness.to(unitReg.inch).magnitude

    xt = (Re - throatRadius)*np.tan(throatTheta)
    absThetaT = abs(throatTheta)

    # thoat arc
    ycTA = throatRadius - (np.sin(np.pi/2 - absThetaT)*designTable["throatArcRad"])
    xcTA = xt - (np.cos(np.pi/2 - absThetaT)*designTable["throatArcRad"])
    arcArr = np.linspace(np.pi/2 - absThetaT, np.pi/2 + np.deg2rad(designTable["convergeAngle"]), circRes)
    throatArc = np.array([nozzle.ContourPoint(xcTA + designTable["throatArcRad"]*np.cos(a), ycTA + designTable["throatArcRad"]*np.sin(a)) for a in arcArr])
    throatArcInner = np.array([nozzle.ContourPoint(xcTA + (designTable["throatArcRad"] - thickness)*np.cos(a), ycTA + (designTable["throatArcRad"] - thickness)*np.sin(a)) for a in arcArr])
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
    convergeArcInner = np.array([nozzle.ContourPoint(xcCA + (designTable["turnArcRad"] + thickness)*np.sin(a), ycCA - (designTable["turnArcRad"] + thickness)*np.cos(a)) for a in arcArr])


    #straight section
    spiketipx = contour[-1].x
    spiketipr = contour[-1].r
    xSS = np.array([x2CL - chamberLength, x2CL - chamberLength, x2CL - chamberLength + thickness, x2CL - chamberLength + thickness])
    rSS = np.array([convergeArc[-1].r, 0, 0, convergeArc[-1].r - thickness])

    straightSection = np.array([nozzle.ContourPoint(xSS[i], rSS[i]) for i in range(len(xSS))])
    
    fullPlugContour = np.concatenate([contour[::-1], throatArc[1:], convergeArc, straightSection, convergeArcInner[::-1], throatArcInner[1:][::-1]], axis=0)
    fullPlugContour = np.insert(fullPlugContour, 0, nozzle.ContourPoint(spiketipx, spiketipr-thickness), axis=0)

    # cooling channels
    baseOuter = np.array([nozzle.ContourPoint(x2CL - chamberLength, convergeArc[-1].r - coolingFloor)])
    baseInner = np.array([nozzle.ContourPoint(x2CL - chamberLength, convergeArc[-1].r - coolingFloor - coolingHeight)])

    arcArr = np.linspace(0, np.deg2rad(designTable["convergeAngle"]), circRes)
    convergeArcOuter = np.array([nozzle.ContourPoint(xcCA + (designTable["turnArcRad"] + coolingFloor)*np.sin(a), ycCA - (designTable["turnArcRad"] + coolingFloor)*np.cos(a)) for a in arcArr])
    convergeArcInner = np.array([nozzle.ContourPoint(xcCA + (designTable["turnArcRad"] + coolingFloor + coolingHeight)*np.sin(a), ycCA - (designTable["turnArcRad"] + coolingFloor + coolingHeight)*np.cos(a)) for a in arcArr])

    arcArr = np.linspace(np.pi/2 + np.deg2rad(designTable["convergeAngle"]), np.pi/2 - absThetaT, circRes)
    throatArcOuter = np.array([nozzle.ContourPoint(xcTA + (designTable["throatArcRad"] - coolingFloor)*np.cos(a), ycTA + (designTable["throatArcRad"] - coolingFloor)*np.sin(a)) for a in arcArr])
    throatArcInner = np.array([nozzle.ContourPoint(xcTA + (designTable["throatArcRad"] - coolingFloor - coolingHeight)*np.cos(a), ycTA + (designTable["throatArcRad"] - coolingFloor - coolingHeight)*np.sin(a)) for a in arcArr])

    nozzleOuter = np.array([])
    nozzleInner = np.array([])
    contourInner = np.array([])

    print(thickness)

    for i in range(len(contour[:-1])):
        angle = np.atan2((contour[i+1].r - contour[i].r),(contour[i+1].x - contour[i].x))
        expandAngle = angle - np.pi/2
        nozzleOuter = np.append(nozzleOuter, nozzle.ContourPoint(contour[i].x + coolingFloor*np.cos(expandAngle), contour[i].r + coolingFloor*np.sin(expandAngle)))
        nozzleInner = np.append(nozzleInner, nozzle.ContourPoint(contour[i].x + (coolingFloor + coolingHeight)*np.cos(expandAngle), contour[i].r + (coolingFloor + coolingHeight)*np.sin(expandAngle)))
        contourInner = np.append(contourInner, nozzle.ContourPoint(contour[i].x + thickness*np.cos(expandAngle), contour[i].r + thickness*np.sin(expandAngle)))

    nozzleOuter = np.append(nozzleOuter, nozzle.ContourPoint(contour[-1].x, contour[-1].r - coolingFloor))
    nozzleInner = np.append(nozzleInner, nozzle.ContourPoint(contour[-1].x, contour[-1].r - coolingFloor - coolingHeight))
    contourInner = np.append(contourInner, nozzle.ContourPoint(contour[-1].x, contour[-1].r - thickness))


    coolOuter = np.concatenate([baseOuter, convergeArcOuter, throatArcOuter, nozzleOuter], axis=0)
    coolInner = np.concatenate([baseInner, convergeArcInner, throatArcInner, nozzleInner], axis=0)

    fullPlugContour = np.append(fullPlugContour, contourInner)

    return fullPlugContour, Q_(abs(xcTA - (x2CL - chamberLength)), unitReg.inch), coolInner, coolOuter

def GenerateDimCowl(throatRadius: Q_, throatTheta: Q_, Re: Q_, straightLength: Q_, chamberOuter: Q_, thickness: Q_, overchoke: Q_, circRes: int = 50):
    designTable = DESIGN.plugDesignTable
    designTable["throatArcRad"] = designTable["throatArcRadFactor"]*Re.to(unitReg.inch).magnitude
    designTable["turnArcRad"] = designTable["turnArcRadFactor"]*Re.to(unitReg.inch).magnitude
    designTable["manifoldDistance"] = designTable["manifoldDistanceFactor"]*Re.to(unitReg.inch).magnitude

    throatRadius = throatRadius.to(unitReg.inch).magnitude
    throatTheta = throatTheta.to(unitReg.radian).magnitude
    Re = Re.to(unitReg.inch).magnitude
    straightLength = straightLength.to(unitReg.inch).magnitude
    chamberOuter = chamberOuter.to(unitReg.inch).magnitude
    thickness = thickness.to(unitReg.inch).magnitude
    coolingFloor = DESIGN.coolingChannelWallDist.to(unitReg.inch).magnitude
    coolingHeight1 = DESIGN.coolingChannelHeightConverge.to(unitReg.inch).magnitude
    coolingHeight2 = DESIGN.coolingChannelHeightChamber.to(unitReg.inch).magnitude
    heightChangeDist = DESIGN.coolingChannelShrinkDist.to(unitReg.inch).magnitude

    absThetaT = abs(throatTheta)
    xt = (Re - throatRadius)*np.tan(throatTheta)
    xc = xt - (np.cos(np.pi/2 - absThetaT)*designTable["throatArcRad"])
    ic(xc, straightLength)

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
    Rmax = chamberOuter + thickness

    xm = designTable["manifoldDistance"]

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

    # chamber wall

    arcArr = np.linspace(-thetaL, np.pi/2 - thetal, circRes)
    lipArc = np.array([nozzle.ContourPoint(xf - rf*np.sin(a), yf - rf*np.cos(a)) for a in arcArr])
    
    # lipStraight = np.array([nozzle.ContourPoint(xf - rf*cl, yf - rf*sl), nozzle.ContourPoint(xc + rc*cl, yc + rc*sl)])

    arcArr = np.linspace(thetal, np.pi/2, circRes)
    convergeArc = np.array([nozzle.ContourPoint(xc + rc*np.cos(a), yc + rc*np.sin(a)) for a in arcArr])

    chamber = np.array([nozzle.ContourPoint(xc - straightLength, yc+rc), nozzle.ContourPoint(xc - straightLength, Rmax)])

    arcArr = np.linspace(np.pi/2, -np.pi/2 + thetaL, circRes)
    manifold = np.array([nozzle.ContourPoint(xm + rm*np.cos(a), ym + rm*np.sin(a)) for a in arcArr])

    final = np.array([nozzle.ContourPoint(0, Re)])

    cowl = np.concatenate([lipArc, convergeArc, chamber, manifold, final], axis=0)

    # cooling

    # arcArr = np.linspace(0, -np.pi/2 + thetaL, circRes)
    arcArr = np.array([-np.pi/2 + thetaL])
    manifold = np.array([nozzle.ContourPoint(xm + (rm - coolingFloor)*np.cos(a), ym + (rm - coolingFloor)*np.sin(a)) for a in arcArr])

    upperManifold = np.array([nozzle.ContourPoint(xm + (rm - coolingFloor - coolingHeight1)*np.cos(a), ym + (rm - coolingFloor - coolingHeight1)*np.sin(a)) for a in arcArr])

    psi = (-thetaL + thetal + np.pi/2)/2 + thetaL
    dist = coolingFloor*np.sqrt(2)
    dist2 = (coolingFloor + coolingHeight1)*np.sqrt(2)
    # ic(psi)
    # arcArr = np.linspace(-thetaL, np.pi/2 - thetal, circRes)
    # lipArc = np.array([nozzle.ContourPoint((xf + dist*(np.cos(psi))) - rf*np.sin(a), (yf + dist*(np.sin(psi))) - rf*np.cos(a)) for a in arcArr])

    lipArc = np.array([nozzle.ContourPoint((xf + dist*(np.cos(psi))) - rf*np.sin(-thetaL), (yf + dist*(np.sin(psi))) - rf*np.cos(-thetaL)),
                       nozzle.ContourPoint((xf + dist*(np.cos(psi))) - rf*np.sin(np.pi/2 - thetal), (yf + dist*(np.sin(psi))) - rf*np.cos(np.pi/2 - thetal))])
    
    lipArcUpper = np.array([nozzle.ContourPoint((xf + dist2*(np.cos(psi))) - rf*np.sin(-thetaL), (yf + dist2*(np.sin(psi))) - rf*np.cos(-thetaL)),
                            nozzle.ContourPoint((xf + dist2*(np.cos(psi))) - rf*np.sin(np.pi/2 - thetal), (yf + dist2*(np.sin(psi))) - rf*np.cos(np.pi/2 - thetal))])

    arcArr = np.linspace(thetal, np.pi/2, circRes)
    convergeArc = np.array([nozzle.ContourPoint(xc + (rc + coolingFloor)*np.cos(a), yc + (rc + coolingFloor)*np.sin(a)) for a in arcArr])

    convergeArcUpper = np.array([nozzle.ContourPoint(xc + (rc + coolingFloor + coolingHeight1)*np.cos(a), yc + (rc + coolingFloor + coolingHeight1)*np.sin(a)) for a in arcArr])

    chamber = np.array([nozzle.ContourPoint(xc - heightChangeDist, yc + (rc + coolingFloor)), nozzle.ContourPoint(xc - straightLength, yc+(rc + coolingFloor))])

    chamberUpper = np.array([nozzle.ContourPoint(xc - heightChangeDist, yc + (rc + coolingFloor + coolingHeight2)), nozzle.ContourPoint(xc - straightLength, yc+(rc + coolingFloor + coolingHeight2))])

    cowlCooling = np.concatenate([manifold, lipArc, convergeArc, chamber], axis=0)
    cowlUpper = np.concatenate([upperManifold, lipArcUpper, convergeArcUpper, chamberUpper], axis=0)

    return cowl, cowlCooling, cowlUpper

def GenerateDimChamber(throatRadius: Q_, throatTheta: Q_, Re: Q_, chamberLength: Q_, chamberOuter: Q_, thickness: Q_, overchoke: Q_, baseRadius: Q_, circRes: int = 50):
    designTable = DESIGN.plugDesignTable
    designTable["throatArcRad"] = designTable["throatArcRadFactor"]*Re.to(unitReg.inch).magnitude
    designTable["turnArcRad"] = designTable["turnArcRadFactor"]*Re.to(unitReg.inch).magnitude
    designTable["manifoldDistance"] = designTable["manifoldDistanceFactor"]*Re.to(unitReg.inch).magnitude

    chamberLength = chamberLength.to(unitReg.inch).magnitude
    baseRadius = baseRadius.to(unitReg.inch).magnitude

    throatRadius = throatRadius.to(unitReg.inch).magnitude
    throatTheta = throatTheta.to(unitReg.radian).magnitude
    Re = Re.to(unitReg.inch).magnitude
    chamberOuter = chamberOuter.to(unitReg.inch).magnitude
    thickness = thickness.to(unitReg.inch).magnitude

    absThetaT = abs(throatTheta)
    xt = (Re - throatRadius)*np.tan(throatTheta)

    origin = np.array([nozzle.ContourPoint(0, Re)])

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

    straightSection = np.array([nozzle.ContourPoint(x2CL - chamberLength, baseRadius)])

    chamber = np.concatenate([throatArc, convergeArc, straightSection], axis=0)

    straightLength = abs(xcTA - (x2CL - chamberLength))

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
    Rmax = chamberOuter + thickness

    xm = designTable["manifoldDistance"]

    Amat = np.array([[0, 0, 1, 0, sL, 0, 0, 0, 0], 
                    [0, 0, 0, 1, -cL, 0, 0, 0 ,0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0], 
                    [0, 0, 0, 0, 0, 1, 1, 0, 0],
                    [0, -cl, 1, 0, -cl, 0, 0, -sl, 0],
                    [1, sl, 0, -1, sl, 0, 0, -cl, 0],
                    [0, 0, 0, 0, 0, 0, -sL, 0, cL],
                    [0, 0, 0, 0, 0, 1, -cL, 0, -sL],
                    [0, 0, cotl, 1, -cscl, 0, 0, 0, 0]])

    bmat = np.array([0, Re, Rinner, Rmax, xcTA, 0, xm, Re, Re-o*cscl])

    yc, rc, xf, yf, rf, ym, rm, l, L = np.linalg.solve(Amat, bmat)
    #! New stuff

    arcArr = np.linspace(-thetaL, np.pi/2 - thetal, circRes)
    lipArc = np.array([nozzle.ContourPoint(xf - rf*np.sin(a), yf - rf*np.cos(a)) for a in arcArr])
    
    arcArr = np.linspace(thetal, np.pi/2, circRes)
    convergeArc = np.array([nozzle.ContourPoint(xcTA + rc*np.cos(a), yc + rc*np.sin(a)) for a in arcArr])

    chamberWall = np.array([nozzle.ContourPoint(xcTA, yc+rc), nozzle.ContourPoint(xcTA - straightLength, yc+rc)])

    wall = np.concatenate([lipArc, convergeArc, chamberWall], axis=0)

    aimpoint = (xcTA, (ycTA + designTable["throatArcRad"] + Rinner)/2)

    return np.concatenate([wall, chamber[::-1], origin], axis=0), aimpoint

def getOverchokeDist(lipRad, throatRad, thetaThroat, areaRatio):
    phi = np.pi/2 + thetaThroat
    return (lipRad - np.sqrt(areaRatio*(lipRad**2 - throatRad**2) + throatRad**2))/np.sin(phi)