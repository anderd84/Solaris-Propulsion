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
    cowlC = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
    chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch), Q_(1.5, unitReg.inch))
    plots.PlotPlug(fig, plugC)
    plots.PlotPlug(fig, cowlC)
    plots.PlotPlug(fig, chamberC)
    # plt.show()

    cooling2 = domain.DomainMC(-8.5, 4, 15, 4, .01)
    cooling2.DefineMaterials(cowlC, np.array([]), chamberC, plugC, 8)

    cooling2.DumpFile("coolmesh.msh")

    # tic = time.perf_counter()
    # cooling = domain.Domain(-8.5, 4, 15, 4, .05, .05)
    # cooling.DefineMaterials(cowlC, np.array([]), chamberC, plugC, fig)
    # toc = time.perf_counter()

    startingpoint = (-8, 2.6) # TODO use real point
    plt.plot([startingpoint[0], aimpoint[0]], [startingpoint[1], aimpoint[1]], 'rx-')
    cooling2.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea, fig)

    # fig2 = plots.CreateNonDimPlot()
    # cooling.ShowMaterialPlot(fig2)

    # cooling2.ShowMaterialPlot(fig)
    cooling2.ShowStatePlot(fig)



    plt.show()

    # print(f"Time single: {toc - tic}")

if __name__ == "__main__":
    main()