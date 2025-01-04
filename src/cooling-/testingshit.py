import matplotlib.pyplot as plt
from Cooling import domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg
from Cooling import cooling2d as cooling_func
from Cooling.material import DomainMaterial 
import numpy as np
from icecream import ic





Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas
print(exhaust.stagTemp)

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

fig = plots.CreateNonDimPlot()
plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch))
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch), Q_(1.5, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
plots.PlotPlug(fig, chamberC)
fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)

coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh2.msh")


                
#TODO Add David's Gamma changing as a function of Temp
coolmesh.ShowStatePlot(fig)
plt.scatter(350 * coolmesh.xstep+ coolmesh.array[0,0].x, -30 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) #top left
plt.scatter(650 * coolmesh.xstep+ coolmesh.array[0,0].x, -30 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) # top right
plt.scatter(350 * coolmesh.xstep+ coolmesh.array[0,0].x, -200 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) #bottome left
plt.scatter(650 * coolmesh.xstep+ coolmesh.array[0,0].x, -200 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) # bottom right
plt.show()