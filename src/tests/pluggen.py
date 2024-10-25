from Nozzle import plug
from fluids import gas
from Nozzle import plots
import matplotlib.pyplot as plt

exhaust = gas.Gas(1.17, 287)

cont, Rt, Tt, field = plug.CreateRaoContour(exhaust, 300, 6200, 6.75, 4.5, 3, 6)
fig = plots.CreateNonDimPlot()
plots.PlotContour(fig, cont, Rt, Tt, 3)
plots.PlotField(fig, field, 3)
plt.show()

