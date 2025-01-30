import numpy as np
import matplotlib.pyplot as plt

from cooling import domain
import general.design as DESIGN
from general.units import Q_, unitReg
from nozzle import plug
from nozzle import plots
from nozzle import analysis

def main():
    Re = Q_(3.2, unitReg.inch)
    exhaust = DESIGN.exhaustGas

    cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
    Rt = outputData["radiusThroat"]
    Tt = outputData["thetaThroat"]
    Re = outputData["radiusLip"]

    fig = plots.CreateNonDimPlot()
    # plots.PlotContour(fig, cont, Rt, Tt, Re)
    # plots.PlotField(fig, field, Re)
    overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)

    plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
    cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
    chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))
    plots.PlotPlug(fig, plugC)
    plots.PlotPlug(fig, cowlC)
    plots.PlotPlug(fig, chamberC)
    fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
    fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
    # plt.show()

    p = Q_(6.75, unitReg.psi)
    rlines, llines, streams = analysis.CalculateComplexField(cont, p, exhaust, 1, Tt, Rt, Re.magnitude, 75, 0, 2)
    fieldGrid = analysis.GridifyComplexField(rlines, llines)

    cooling2 = domain.DomainMC(-7.3, 4.1, 7.9, 3, .05)
    cooling2.DefineMaterials(cowlC, chamberC, plugC, 10)


    # tic = time.perf_counter()
    # cooling = domain.Domain(-8.5, 4, 15, 4, .05, .05)
    # cooling.DefineMaterials(cowlC, np.array([]), chamberC, plugC, fig)
    # toc = time.perf_counter()

    startingpoint = (-6.75, 2.6) # TODO use real point
    plt.plot([startingpoint[0], aimpoint[0]], [startingpoint[1], aimpoint[1]], 'rx-')

    cooling2.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea)
    # cooling2.AssignExternalTemps(fieldGrid, exhaust, DESIGN.chokeArea)
    # cooling2.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(400, unitReg.psi))
    # cooling2.AssignCoolantFlow(domain.CoolingChannel(plugCoolU, plugCoolL), True, Q_(400, unitReg.psi))

    cooling2.DumpFile("coolmesh")

    # fig2 = plots.CreateNonDimPlot()
    # cooling2.ShowMaterialPlot(fig, False)

    # cooling2.ShowMaterialPlot(fig)
    cooling2.NodePlot(fig, "temperature")

    fig2 = plots.CreateNonDimPlot()
    analysis.PlotFieldData(fig2, fieldGrid, 1, 1)

    plt.show()

    # print(f"Time single: {toc - tic}")

if __name__ == "__main__":
    main()