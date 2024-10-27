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
plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re)
# plt.show()

cooling = domain.Domain(-.5, 3.5, 6, 3.5, 500, 500)
cooling.DefineMaterials(np.array([]), np.array([]), np.array([]), cont, fig)
cooling.ShowMaterialPlot(fig)
plt.show()
