from Nozzle import plug
from fluids import gas
from Nozzle import plots
from Nozzle import analysis
import matplotlib.pyplot as plt
from icecream import ic
import numpy as np
np.product = np.prod
import General.design as DESIGN
from General.units import Q_, unitReg
import matrix_viewer as mv

# import pint

# ureg = pint.UnitRegistry()
# ureg.default_system = 'US'
# Q_ = ureg.Quantity

# exhaustt = gas.Gas(1.17, 287)
# exhaustt.Rgas = Q_(68.004, ureg.foot * ureg.force_pound / ureg.degR / ureg.pound)
# exhaustt.stagPress = Q_(300, ureg.psi)
# exhaustt.stagTemp = Q_(6200, ureg.degR)
# ic(exhaustt.getChokedArea(Q_(8, ureg.pound/ureg.second)).to(ureg.inch**2))


Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas
ic(exhaust.stagTemp.to(unitReg.degR))
ic(exhaust.stagPress.to(unitReg.psi))
ic(exhaust.Rgas.to(unitReg.joule/unitReg.kg/unitReg.kelvin))
ic(DESIGN.totalmdot.to(unitReg.kg/unitReg.second))

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]
ic(outputData["areaRatio"])
ic(np.rad2deg(Tt))
ic(Re)
phi = np.pi/2 + Tt
Astar = np.pi/np.sin(phi) * (Re**2 - Rt**2)
ic(Astar)
ic(DESIGN.chokeArea)

Cf = outputData["Cf"]
thrust = Astar * DESIGN.chamberPressure * Cf
ic(thrust)
ic(Cf)



fig = plots.CreateNonDimPlot()
plots.PlotContour(fig, cont, Rt, Tt, Re)
plt.plot([cont[-1].x], [cont[-1].r], 'ro')
# plots.PlotField(fig, field, Re)
# plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
# cowlC = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
# plots.PlotPlug(fig, plugC)
# plots.PlotPlug(fig, cowlC)
rlines, llines, streams = analysis.CalculateComplexField(cont, Q_(.01, unitReg.psi), exhaust, 1, Tt, Rt, Re.magnitude, 100, 0, 1)
istream = streams[0]
fig.axes[0].plot([p.x for p in istream], [p.r for p in istream], '--b', linewidth=1.5)
ostream = streams[1]
fig.axes[0].plot([p.x for p in ostream], [p.r for p in ostream], '--b', linewidth=1.5)
fieldGrid = analysis.GridifyComplexField(rlines, llines)
analysis.PlotFieldData(fig, fieldGrid)
analysis.PlotCharacteristicLines(fig, np.concatenate((rlines, llines), axis=0))

x = np.array([[p.x for p in row] for row in rlines])
r = np.array([[p.r for p in row] for row in rlines])
term = np.array([[p.terminate for p in row] for row in rlines])

# mv.view(x)
# mv.view(r)
# mv.view(term)

# mv.show()

analysis.CalculateThrust(exhaust, Q_(1, unitReg.psi), Tt, Rt, Re, istream, cont[-1].r)
plt.show()

# plots.WriteContourTXT(plugC, "plug.txt")