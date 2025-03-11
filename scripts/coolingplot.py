from alive_progress import alive_bar
from cooling import material
from cooling.material import CoolantType, DomainMaterial, MaterialType
import matplotlib.pyplot as plt
import numpy as np

from cooling import domain
from nozzle import plots
import general.design as DESIGN
from nozzle import plug
from nozzle import analysis
from general.units import Q_, unitReg

Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

# fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re)
overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)

plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))

fig2 = plots.CreateNonDimPlot()

plots.PlotPlug(fig2, plugC)
plots.PlotPlug(fig2, cowlC)
plots.PlotPlug(fig2, chamberC)
fig2.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
fig2.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
fig2.axes[0].plot([p.x for p in plugCoolL], [p.r for p in plugCoolL], '-k', linewidth=1)
fig2.axes[0].plot([p.x for p in plugCoolU], [p.r for p in plugCoolU], '-k', linewidth=1)

mesh = domain.DomainMC.LoadFile("save")
# mesh = domain.DomainMC.LoadFile("highestmesh2")

# for i in range(mesh.vpoints):
#     for j in range(mesh.hpoints):
#         mesh.array[i,j].temperature -= mesh2.array[i,j].temperature

# shitOnes = set()
# for i in range(mesh.vpoints):
#     for j in range(mesh.hpoints):
#         if mesh.array[i,j].temperature.m < 0:
#             shitOnes.add((i,j))
#             plt.plot(mesh.array[i,j].x, mesh.array[i,j].r, 'go')
#         if mesh.array[i,j].temperature.m > 1e5:
#             shitOnes.add((i,j))
#             plt.plot(mesh.array[i,j].x, mesh.array[i,j].r, 'bo')

# print(shitOnes)

mesh.NodePlot(fig2, "temperature", [DomainMaterial.FREE, DomainMaterial.CHAMBER, DomainMaterial.EXHAUST])
mesh.NodePlot(plots.CreateNonDimPlot(), "temperature", [DomainMaterial.FREE, DomainMaterial.CHAMBER, DomainMaterial.EXHAUST, DomainMaterial.COWL, DomainMaterial.PLUG])
mesh.NodePlot(plots.CreateNonDimPlot(), "pressure", [DomainMaterial.FREE, DomainMaterial.CHAMBER, DomainMaterial.EXHAUST, DomainMaterial.COWL, DomainMaterial.PLUG])
mesh.NodePlot(plots.CreateNonDimPlot(), "material", [DomainMaterial.FREE, DomainMaterial.CHAMBER, DomainMaterial.EXHAUST])

# mesh.RelationPlot(fig2)
# mesh.RelationPlot(fig2)

plt.show()
