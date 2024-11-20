from General import units
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


overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)
ic(overchoke)


fig = plots.CreateNonDimPlot()
plots.PlotContour(fig, cont, Rt, Tt, Re)
plt.plot([cont[-1].x], [cont[-1].r], 'ro')
# plots.PlotField(fig, field, Re)
plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCool, asv = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
chamberC = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))
ic(chamberC)
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
# plots.PlotPlug(fig, chamberC)
plt.plot([p.x for p in cowlCool], [p.r for p in cowlCool], '-b', linewidth=1.5)
plt.plot([p.x for p in asv], [p.r for p in asv], '-r', linewidth=1.5)
# mat, stream = analysis.CalculateComplexField(cont, Q_(1, unitReg.psi), exhaust, 1, Tt, Re.magnitude, 25, 0, 1)
# fig.axes[0].plot([p.x for p in stream], [p.r for p in stream], '--b', linewidth=1.5)
# fieldGrid = analysis.GridifyComplexField(mat, np.array([]))
# analysis.PlotFieldData(fig, fieldGrid)
# analysis.PlotCharacteristicLines(fig, mat)
# analysis.CalculateThrust(exhaust, Q_(1, unitReg.psi), Tt, Rt, Re, fieldGrid, cont[-1])
plt.show()

# plots.WriteContourTXT(plugC, "plug.txt")