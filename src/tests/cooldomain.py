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
plots.PlotField(fig, field, Re)
plt.show()

cooling = domain.Domain(0, 4, 10, 4, 10, 10)
cooling.