from Nozzle import plug
from fluids import gas
from Nozzle import plots
import matplotlib.pyplot as plt
from icecream import ic

exhaust = gas.Gas(1.17, 287)

cont, field, outputData = plug.CreateRaoContour(exhaust, 300, 6200, 6.75, 45, 3, 8)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Cf = outputData["Cf"]

ic(Cf)
ic(outputData["thetaLip"])

fig = plots.CreateNonDimPlot()
plots.PlotContour(fig, cont, Rt, Tt, 3)
plots.PlotField(fig, field, 3)
plt.show()
