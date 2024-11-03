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

Re = Q_(3, unitReg.inch)
exhaust = DESIGN.exhaustGas

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]

fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
plugC = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC = plug.GenerateDimCowl(Rt, Tt, Re, Q_(8, unitReg.inch), DESIGN.maxRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
# plt.show()

tic2 = time.perf_counter()
cooling2 = domain.DomainMC(-8.5, 4, 15, 4, .05)
cooling2.DefineMaterials(cowlC, np.array([]), np.array([]), plugC, fig)
toc2 = time.perf_counter()

# tic = time.perf_counter()
# cooling = domain.Domain(-8.5, 4, 15, 4, .05, .05)
# cooling.DefineMaterials(cowlC, np.array([]), np.array([]), plugC, fig)
# toc = time.perf_counter()

cooling2.ShowMaterialPlot(fig)

# fig2 = plots.CreateNonDimPlot()
# cooling.ShowMaterialPlot(fig2)

cooling2.DumpFile("coolmesh.pkl")

plt.show()

# print(f"Time single: {toc - tic}")
print(f"Time multi: {toc2 - tic2}")
