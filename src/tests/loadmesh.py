import matplotlib.pyplot as plt

from Cooling import domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg

Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas
print(exhaust.stagTemp)

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re)
plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch))
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch), Q_(1.5, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
plots.PlotPlug(fig, chamberC)
# fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)


coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh.msh")

startingpoint = (-5.75, 2.6) # TODO use real point
# plt.plot([startingpoint[0], aimpoint[0]], [startingpoint[1], aimpoint[1]], 'rx-')
# coolmesh.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea, fig)

coolmesh.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(400, unitReg.psi))

# fig2 = plots.CreateNonDimPlot()
# coolmesh.ShowMaterialPlot(fig)

# cooling2.ShowMaterialPlot(fig)
# cooling2.ShowStatePlot(fig)
coolmesh.ShowBorderPlot(fig)


plt.show()