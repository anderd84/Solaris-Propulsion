import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.optimize import fsolve
import gas
import nozzle

def CalculateSimpleField(contour, PambPc, PbPc, gamma, Mt, Tt, steps = 100, reflections = 3):
    phit = np.pi/2 + Tt
    Me = np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2])
    machs = np.linspace(Mt, Me, steps)
    temps = 1/(1 + gamma[2]*machs**2)
    press = temps**(gamma[5])
    charLineAngles = Tt + (gas.PrandtlMeyerFunction(machs, gamma) - gas.PrandtlMeyerFunction(Mt, gamma)) - gas.MachAngle(machs)
    thetaEnd = Tt + (gas.PrandtlMeyerFunction(Me, gamma) - gas.PrandtlMeyerFunction(Mt, gamma))

    nu1 = gas.PrandtlMeyerFunction(Mt, gamma)
    nu2 = gas.PrandtlMeyerFunction(Me, gamma)

    print(f"Me: {Me}, Mt: {Mt}")
    print(f"ThetaEnd: {np.rad2deg(thetaEnd)}, Tt: {np.rad2deg(Tt)}")

    endpoints, bounds = CalculateEndpoints(contour, charLineAngles)
    # gas.weakObliqueShock(2, np.deg2rad(10), gamma)
    temps = temps[bounds[0]:bounds[1]]
    machs = machs[bounds[0]:bounds[1]]
    # dA = CalculateDiffArea(endpoints)


    plt.figure()
    plt.plot([p.x for p in endpoints], [p.r for p in endpoints], '.r')
    plt.plot([p.x for p in contour], [p.r for p in contour], '--k', linewidth=0.5)
    # xgrid, ygrid, chord = SimpleCharacteristicsMesh(endpoints, machs)
    # plt.contourf(xgrid, ygrid, chord, levels=100, cmap='jet')

    ColorLine(endpoints[:], temps[:], 'jet', plt.gca())
    plt.plot([1, 1 + np.cos(thetaEnd)], [1, 1 + np.sin(thetaEnd)], '-r', linewidth=2)
    plt.plot([0, np.sin(Tt)], [1, 1 - np.cos(Tt)], '-b', linewidth=2)

    plt.plot([0, 10*np.cos(Tt + nu1)], [1, 1+10*np.sin(Tt + nu1)], '--r', linewidth=1)
    plt.plot([0, 10*np.cos(Tt + nu2 - nu1)], [1, 1+10*np.sin(Tt + nu2 - nu1)], '--g', linewidth=1)

    for ep in endpoints:
        plt.plot([0, ep.x], [1, ep.r], '--k', linewidth=0.1)
    


    plt.xlim(contour[0].x, contour[-1].x)
    plt.ylim(0, 2)
    plt.gca().set_aspect('equal')
    plt.show()


def CalculateEndpoints(contour, angles):
    rayLength = np.sqrt((0 - contour[-1].x)**2 + (1 - contour[-1].r)**2)
    
    endpoints = np.empty(len(angles), dtype=nozzle.ContourPoint)

    lowerBound = 0
    for i, angle in enumerate(angles):
        for j in range(len(contour) - 1):
            a, b, c, d = contour[j].x, contour[j].r, contour[j+1].x, contour[j+1].r
            rayEndx = 0 + rayLength * np.cos(angle)
            rayEndr = 1 + rayLength * np.sin(angle)
            x = ((d - (b - d)/(a - c)*c) - 1)/((rayEndr - 1)/(rayEndx) - (b - d)/(a - c))
            if a <= x <= c:
                endpoints[i] = nozzle.ContourPoint(x, (rayEndr - 1)/(rayEndx)*x + 1)
                break
            if j == len(contour) - 2: 
                if endpoints[i - 1] is not None:
                    return endpoints[lowerBound:i], (lowerBound, i)
                else:
                    lowerBound = i + 1
    return endpoints, (lowerBound, len(endpoints))

def ColorLine(points, colorFunc, colormap, axes):
    x, y = [p.x for p in points], [p.r for p in points]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=colormap, norm=plt.Normalize(min(colorFunc), max(colorFunc)), array=colorFunc)
    line = axes.add_collection(lc)
    plt.colorbar(line, ax=axes)
    return line

def SimpleCharacteristicsMesh(contour, chordFunc):
    xgrid = np.array([p.x for p in contour])
    xgrid = np.array([np.zeros(len(contour)), xgrid])

    rgrid = np.array([p.r for p in contour])
    rgrid = np.array([np.ones(len(contour)), rgrid])

    chord = np.array([chordFunc[:len(contour)], chordFunc[:len(contour)]])

    return xgrid, rgrid, chord

def CalculateDiffArea(contour):
    dA = np.zeros(len(contour) - 1)
    for i in range(len(contour) - 1):
        a, b, c, d = contour[i].x, contour[i].r, contour[i+1].x, contour[i+1].r
        theta = np.arctan((d - b)/(c - a))
        dA[i] = abs(np.pi/np.sin(theta)*(b**2 - d**2))
    return dA
