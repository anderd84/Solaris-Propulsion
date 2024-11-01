from Nozzle import plug
from fluids import gas
from Nozzle import plots
import matplotlib.pyplot as plt
from icecream import ic
import numpy as np
# import pint

# ureg = pint.UnitRegistry()
# ureg.default_system = 'US'
# Q_ = ureg.Quantity

# exhaustt = gas.Gas(1.17, 287)
# exhaustt.Rgas = Q_(68.004, ureg.foot * ureg.force_pound / ureg.degR / ureg.pound)
# exhaustt.stagPress = Q_(300, ureg.psi)
# exhaustt.stagTemp = Q_(6200, ureg.degR)
# ic(exhaustt.getChokedArea(Q_(8, ureg.pound/ureg.second)).to(ureg.inch**2))


Re = 3
exhaust = gas.Gas(1.17, 287)

cont, field, outputData = plug.CreateRaoContour(exhaust, 300, 6200, 6.75, 15, Re, 5)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
ic(outputData["areaRatio"])
ic(np.rad2deg(Tt))
phi = np.pi/2 + Tt
Astar = np.pi/np.sin(phi) * (Re**2 - Rt**2)
ic(Astar)
Cf = outputData["Cf"]
thrust = Astar * 300 * Cf
ic(thrust)
ic(Cf)
ic(outputData["thetaLip"])

fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re)
plugC = plug.GenerateDimPlug(cont, Rt, Tt, Re, 8, 2)
cowlC = plug.GenerateDimCowl(cont, Rt, Tt, Re, 8, 3.5, .25)
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
mat, stream = analysis.CalculateComplexField(cont, 14/300, 15/300, exhaust, 1, Tt, Re.magnitude, 25, 0, 3)
fig.axes[0].plot([p.x for p in stream], [p.r for p in stream], '--b', linewidth=1.5)
fieldGrid = analysis.GridifyComplexField(mat, np.array([]))
analysis.PlotFieldData(fig, fieldGrid)
analysis.PlotCharacteristicLines(fig, mat)
plt.show()
