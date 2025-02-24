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
print(exhaust.stagTemp)

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

throatHyroD = 2*(Re - Rt)/np.cos(Tt)
# fig = plots.CreateNonDimPlot()
# plots.PlotContour(fig, cont, Rt, Tt, Re)
# plots.PlotField(fig, field, Re)
overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)

plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))
# plots.PlotPlug(fig, plugC)
# plots.PlotPlug(fig, cowlC)
# plots.PlotPlug(fig, chamberC)
# fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in plugCoolL], [p.r for p in plugCoolL], '-k', linewidth=1)
# fig.axes[0].plot([p.x for p in plugCoolU], [p.r for p in plugCoolU], '-k', linewidth=1)

# plt.show()


# # mmapmesh = domain.DomainMMAP(coolmesh)
# coolmesh.ShowStatePlot(fig, "temperature")

startingpoint = (-6.75, 2.6) # TODO use real point


# highmesh = domain.DomainMC.LoadFile("highmesh")

highmesh = domain.DomainMC(-7.3, 4.1, 13.5, 3.75, .01)
p = Q_(6.75, unitReg.psi)
rlines, llines, streams = analysis.CalculateComplexField(cont, p, exhaust, 1, Tt, Rt, Re.magnitude, 75, 5, 2)
fieldGrid = analysis.GridifyComplexField(rlines, llines)
# analysis.PlotFieldData(fig, fieldGrid, 1, 1)

highmesh.DefineMaterials(cowlC, chamberC, plugC, 10)
highmesh.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea)
highmesh.AssignExternalTemps(fieldGrid, cont, exhaust, DESIGN.chokeArea, throatHyroD)

# coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("save")
# highmesh.ApplyStateMap(coolmesh, {"temperature", "pressure"})

highmesh.DumpFile("highmesh")
fig2 = plots.CreateNonDimPlot()

outerloop = highmesh.NewCoolantLoop(Q_(.03, 'inch'), 120, DESIGN.Fuel_Total, CoolantType.RP1)
highmesh.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(360, unitReg.psi), outerloop)
innerloop = highmesh.NewCoolantLoop(Q_(.03, 'inch'), 90, Q_(2, unitReg.pound/unitReg.sec), CoolantType.RP1)
highmesh.AssignCoolantFlow(domain.CoolingChannel(plugCoolU, plugCoolL), True, Q_(100, unitReg.psi), innerloop)

# print(highmesh.array[0,0])
highmesh.GuessChannelState(outerloop, Q_(650, unitReg.degR))
highmesh.GuessChannelState(innerloop, Q_(650, unitReg.degR))
print(highmesh.coolingLoops)
highmesh.DumpFile("highmesh2")
plots.PlotPlug(fig2, plugC)
plots.PlotPlug(fig2, cowlC)
plots.PlotPlug(fig2, chamberC)
fig2.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
fig2.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
fig2.axes[0].plot([p.x for p in plugCoolL], [p.r for p in plugCoolL], '-k', linewidth=1)
fig2.axes[0].plot([p.x for p in plugCoolU], [p.r for p in plugCoolU], '-k', linewidth=1)

calcPoints = set()
blacklist = set()
with alive_bar(highmesh.vpoints*highmesh.hpoints, title="Finding calculation points") as bar:
    for row in range(highmesh.vpoints):
        for col in range(highmesh.hpoints):
            if highmesh.array[row, col].material not in material.MaterialType.STATIC_TEMP:
                calcPoints.add((row, col))
            if highmesh.array[row, col].material == material.DomainMaterial.COOLANT_BULK and highmesh.array[row, col].border:
                blacklist.add(highmesh.array[row, col].previousFlow)
            bar()

for pair in blacklist:
    calcPoints.remove(pair)

plotx = [highmesh.array[pnt].x for pnt in calcPoints]
plotr = [highmesh.array[pnt].r for pnt in calcPoints]

# plt.plot(plotx, plotr, 'go')

highmesh.NodePlot(fig2, "temperature", [DomainMaterial.CHAMBER, DomainMaterial.FREE])#, DomainMaterial.PLUG, DomainMaterial.COWL, DomainMaterial.COOLANT_INLET, DomainMaterial.COOLANT_OUTLET])
# highmesh.RelationPlot(fig2)
# highmesh.ShowCellPlot(fig2)

# cooling2.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(400, unitReg.psi))
# cooling2.AssignCoolantFlow(domain.CoolingChannel(plugCoolU, plugCoolL), True, Q_(400, unitReg.psi))
# mmapmesh.plotPressDrop(fig)
# mv.view(np.asarray(mmapmesh.pressure.magnitude))
# mv.show()

# plt.plot([startingpoint[0], aimpoint[0]], [startingpoint[1], aimpoint[1]], 'rx-')
# coolmesh.AssignChamberTemps(chamberC, exhaust, startingpoint, aimpoint, DESIGN.chamberInternalRadius, DESIGN.plugBaseRadius, DESIGN.chokeArea, fig)

# coolmesh.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(400, unitReg.psi))
# coolmesh.AssignCoolantFlow(domain.CoolingChannel(plugCoolU, plugCoolL), True, Q_(400, unitReg.psi))

# coolmesh.DumpFile("coolmesh.msh")

# coolmesh.ShowStatePlot(fig, "temperature")
# coolmesh.ShowMaterialPlot(fig)

# coolmesh.ShowBorderPlot(fig)


plt.show()