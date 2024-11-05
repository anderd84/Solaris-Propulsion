import matplotlib.pyplot as plt

from Cooling import domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg
import sys
from pympler.asizeof import asizeof

Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re)
plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch), Q_(1.5, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
plots.PlotPlug(fig, chamberC, '-r')


coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh.pkl")

startingpoint = (-8, 2.6) # TODO use real point
plt.plot([startingpoint[0], aimpoint[0]], [startingpoint[1], aimpoint[1]], 'rx-')
coolmesh.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.chokeArea, fig)

# fig2 = plots.CreateNonDimPlot()
# cooling.ShowMaterialPlot(fig2)

# coolmesh.ShowMaterialPlot(fig)
coolmesh.ShowStatePlot(fig)

# coolmesh.ShowMaterialPlot(fig)


plt.show()