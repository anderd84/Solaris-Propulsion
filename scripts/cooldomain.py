from cooling import domain
from nozzle import plug
from nozzle import plots
import general.design as DESIGN
from general.units import Q_, unitReg

import numpy as np
import matplotlib.pyplot as plt

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

cooling = domain.Domain(-8.5, 4, 15, 4, .1, .1)
cooling.DefineMaterials(cowlC, np.array([]), np.array([]), plugC, fig)
cooling.ShowMaterialPlot(fig)
plt.show()
