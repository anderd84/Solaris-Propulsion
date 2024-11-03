import matplotlib.pyplot as plt

from Cooling import domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg
import sys
from pympler.asizeof import asizeof

Re = Q_(3, unitReg.inch)
exhaust = DESIGN.exhaustGas

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]

fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
plugC = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC = plug.GenerateDimCowl(Rt, Tt, Re, Q_(8, unitReg.inch), DESIGN.maxRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)


coolmesh = domain.DomainMC.LoadFile("coolmesh.pkl")

print(asizeof(coolmesh))

# coolmesh.ShowMaterialPlot(fig)
print(f"{coolmesh.array.nbytes / 1000000} MB" )
coolmesh.DumpFile("coolmesh.pkl")

plt.show()