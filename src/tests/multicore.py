from Cooling import domain
from Cooling import material
from fluids import gas
from Nozzle import plug
from Nozzle import plots
import General.design as DESIGN
from General.units import Q_, unitReg
import time

import numpy as np
import matplotlib.pyplot as plt
from icecream import ic

import pickle
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
    plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
    cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch))
    chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch), Q_(1.5, unitReg.inch))
    plots.PlotPlug(fig, plugC)
    plots.PlotPlug(fig, cowlC)
    plots.PlotPlug(fig, chamberC)
    fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
    fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
    # plt.show()

    cooling2 = domain.DomainMC(-6, 4.5, 7, 5, .01)
    cooling2.DefineMaterials(cowlC, np.array([]), chamberC, plugC, 17)


    # tic = time.perf_counter()
    # cooling = domain.Domain(-8.5, 4, 15, 4, .05, .05)
    # cooling.DefineMaterials(cowlC, np.array([]), chamberC, plugC, fig)
    # toc = time.perf_counter()

    startingpoint = (-5.75, 2.6) # TODO use real point
    plt.plot([startingpoint[0], aimpoint[0]], [startingpoint[1], aimpoint[1]], 'rx-')
    cooling2.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea, fig)

    cooling2.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(400, unitReg.psi))

    cooling2.DumpFile("coolmesh.msh")

    # fig2 = plots.CreateNonDimPlot()
    cooling2.ShowMaterialPlot(fig)

    # cooling2.ShowMaterialPlot(fig)
    # cooling2.ShowStatePlot(fig)



    plt.show()

    # print(f"Time single: {toc - tic}")

if __name__ == "__main__":
    main()