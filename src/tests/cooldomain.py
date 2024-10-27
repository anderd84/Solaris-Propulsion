from Cooling import domain
from Cooling import material
from fluids import gas
from Nozzle import plug
from Nozzle import plots

import numpy as np
import matplotlib.pyplot as plt
from icecream import ic


Re = 3
exhaust = gas.Gas(1.17, 287)

cont, field, outputData = plug.CreateRaoContour(exhaust, 300, 6200, 6.75, 15, Re, 5)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]

fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
plugC = plug.GenerateDimPlug(cont, Rt, Tt, Re, 8, 1)
cowlC = plug.GenerateDimCowl(cont, Rt, Tt, Re, 8, 3.5, .25)
# plots.PlotPlug(fig, plugC)
# plots.PlotField(fig, field, Re)
# plt.show()

cooling = domain.Domain(-8.5, 4, 15, 4, .025, .025)
cooling.DefineMaterials(cowlC, np.array([]), np.array([]), plugC, fig)
cooling.ShowMaterialPlot(fig)
plt.show()
