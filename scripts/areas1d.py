import dataclasses
import numpy as np
import matplotlib.pyplot as plt
import pint
from alive_progress import alive_it
import pandas as pd
from scipy.optimize import fsolve

from cooling import domain
from cooling.material import DomainMaterial
from fluids import gas
import general.design as DESIGN
from general.units import Q_, unitReg
from nozzle import plug
from nozzle import plots
from nozzle import analysis

def main():
    Re = Q_(3.2, unitReg.inch)
    exhaust = DESIGN.exhaustGas
    fig = plots.CreateNonDimPlot()

    cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
    Rt = outputData["radiusThroat"]
    Tt = outputData["thetaThroat"]
    Re = outputData["radiusLip"]

    overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)

    plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
    cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
    chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))
    plots.PlotPlug(fig, plugC)
    plots.PlotPlug(fig, cowlC)
    plots.PlotPlug(fig, chamberC)
    fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
    fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)


    cooling2 = domain.DomainMC(-7.3, 4.1, 7.9, 3, .05)
    cooling2.DefineMaterials(cowlC, chamberC, plugC, 10)

    startingpoint = (-6.75, 2.6) # TODO use real point

    cooling2.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea)

    edgeData = []

    for i in alive_it(range(len(chamberC) - 1)):
        dist = np.sqrt((chamberC[i].x - chamberC[i+1].x)**2 + (chamberC[i].r - chamberC[i+1].r)**2)

        steps = max(int(dist/min(cooling2.xstep, cooling2.rstep) * 1.5), 5)

        xpts = np.linspace(chamberC[i].x, chamberC[i+1].x, steps)[:-1]
        rpts = np.linspace(chamberC[i].r, chamberC[i+1].r, steps)[:-1]

        for j in range(steps - 1):
            if xpts[j] > cooling2.x0 + cooling2.width or xpts[j] < cooling2.x0:
                break
            if rpts[j] < cooling2.r0 - cooling2.height or rpts[j] > cooling2.r0:
                break

            row, col = cooling2.CoordsToCell(xpts[j], rpts[j])
            if cooling2.array[row, col].material != DomainMaterial.CHAMBER:
                offsets = [(0, 1), (0, -1), (1, 0), (-1, 0)]
                for offset in offsets:
                    if cooling2.array[row + offset[0], col + offset[1]].material == DomainMaterial.CHAMBER:
                        row = row + offset[0]
                        col = col + offset[1]
                        break
                else:
                    continue
            wallData = dataclasses.replace(cooling2.array[row, col])
            wallData.x = xpts[j]
            wallData.r = rpts[j]
            edgeData.append(wallData)


    # print(edgePoints)
    # print(edgePoints)

    x = np.array([point.x for point in edgeData])
    r = np.array([point.r for point in edgeData])
    temp = np.array([point.temperature.m for point in edgeData])
    pressure = np.array([point.pressure.m for point in edgeData])
    area = np.array([point.area.m for point in edgeData])
    flowHeight = np.array([point.flowHeight.m for point in edgeData])
    hydroD = np.array([point.hydraulicDiameter.m for point in edgeData])
    velocity = np.array([point.velocity.m for point in edgeData])
    astar = DESIGN.chokeArea.m
    AAstar = area/astar
    mach = np.array([])
    for aa in AAstar:
        machp = fsolve(lambda M: gas.Isentropic1DExpansion(M, exhaust.gammaTyp) - aa, .25)[0]
        machp = min(1, machp)
        mach = np.append(mach, machp)

    df = pd.DataFrame({
        "x": x - min(x),
        "r": r,
        "temp": temp,
        "press": pressure,
        "area": area,
        "flowHeight": flowHeight,
        "hydroD": hydroD,
        "vel": velocity,
        "mach": mach
    })
    df = df.sort_values("x", ascending=True)
    df.to_excel("edgeData.xlsx")
    print(df)
    plt.plot(x, r, 'ro')
    cooling2.ShowCellPlot(fig)

    plt.show()

if __name__ == "__main__":
    main()