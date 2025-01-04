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

ic(gas.MachToVelocity(2.916, exhaust).to(unitReg.foot/unitReg.second))

fig = plots.CreateNonDimPlot()
plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re.magnitude)
# plt.plot([cont[-1].x], [cont[-1].r], 'ro')
# plots.PlotField(fig, field, Re)
plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
# chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))
# plots.PlotPlug(fig, plugC)
# plots.PlotPlug(fig, cowlC)
# plots.PlotPlug(fig, chamberC)
# fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in plugCoolL], [p.r for p in plugCoolL], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in plugCoolU], [p.r for p in plugCoolU], '-k', linewidth=1)
# print(units.PRESCOTT_PRESSURE)
p = Q_(5, unitReg.psi)
# p = units.PRESCOTT_PRESSURE
ic(p.to(unitReg.psi))
rlines, llines, streams = analysis.CalculateComplexField(cont, p, exhaust, 1, Tt, Rt, Re.magnitude, 25, 0, 2)
istream = streams[0]
fig.axes[0].plot([p.x for p in istream], [p.r for p in istream], '--b', linewidth=1.5)
ostream = streams[1]
fig.axes[0].plot([p.x for p in ostream], [p.r for p in ostream], '--b', linewidth=1.5)
fieldGrid = analysis.GridifyComplexField(rlines, llines)

analysis.PlotFieldData(fig, fieldGrid, 20, 20)
# analysis.PlotCharacteristicLines(fig, np.concatenate((rlines, llines), axis=0))

T = (analysis.CalculateThrust(exhaust, p, Tt, Rt, Re, istream, cont[-1].r).to(unitReg.pound_force))
ic(T)
ic((T/DESIGN.totalmdot/Q_(32.2, unitReg.ft/unitReg.s/unitReg.s)).to(unitReg.s))


g0 = Q_(32.2, unitReg.feet/unitReg.second/unitReg.second)
janusVacIsp = Q_(304, unitReg.second)
janusmdot = Q_(7.5, unitReg.pound/unitReg.seconds)
janusVacThrust = janusVacIsp * janusmdot * g0
janusAe = np.pi * (Q_(6.25, unitReg.inch)/2)**2

janusThrust = janusVacThrust - p*janusAe
ic(janusThrust.to(unitReg.pound_force))
plt.show()

# plots.WriteContourTXT(cont[::-1], "cont.txt")
