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

coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh.msh")


total_change = 0
for iterate in range(1):
    ic(iterate)
#    fig.axes[0].clear()
#    coolmesh.ShowStatePlot(fig)
#    fig.canvas.draw()
#    fig.canvas.flush_events()
    #for i in range(coolmesh.vpoints):
    for i in range(30, 200):
        #for j in range(coolmesh.hpoints):
        for j in range(350, 650):
            #*Finding all options for barrier
            if i==122 and j == 578:
                pass
            if not(coolmesh.array[i,j].material == DomainMaterial.CHAMBER or coolmesh.array[i,j].material == DomainMaterial.FREE):  # Select only walls and coolant
                
                if (coolmesh.array[i,j].material == DomainMaterial.COOLANT_WALL or coolmesh.array[i,j].material == DomainMaterial.COOLANT_BULK):    # Select coolant
                        T_new = coolant(coolmesh,i,j)
                        coolmesh.array[i,j].temperature = Q_(T_new.magnitude, unitReg.degR)
                        if T_new.magnitude > 600:
                            pass
                        continue    # Move to next iteration
                if not(coolmesh.array[i,j].border): # Select non-border coolant nodes
                    C_left, C_upper, C_bottom, C_right, T_left, T_upper, T_bottom, T_right = getcorecond(coolmesh,i,j)
                else:   # Wall nodes
   
                    C_left, T_left, C_right, T_right = horizontalcond(coolmesh,i,j)
                    C_upper, T_upper, C_bottom, T_bottom = verticalcond(coolmesh,i,j)
                T_left, T_right, T_upper, T_bottom = Q_([T_left.magnitude, T_right.magnitude, T_upper.magnitude, T_bottom.magnitude], unitReg.degR)#
                Num = (C_left * T_left + C_upper * T_upper + C_bottom * T_bottom + C_right*  T_right)
                Denom = (C_left + C_upper + C_bottom + C_right)
                new_temp = (Num/Denom).to(unitReg.degR)
                current_temp = coolmesh.array[i,j].temperature.to(unitReg.degR)
                current_temp = Q_(current_temp.magnitude, unitReg.degR)
                total_change = total_change + np.abs(current_temp - new_temp)
                check2 = new_temp
                coolmesh.array[i,j].temperature = new_temp
                if new_temp.magnitude > 600:
                    pass
            



                
#TODO Add David's Gamma changing as a function of Temp
#coolmesh.ShowStatePlot(fig)
coolmesh.ShowMaterialPlot(fig)
plt.scatter(350 * coolmesh.xstep+ coolmesh.array[0,0].x, -30 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) #top left
plt.scatter(650 * coolmesh.xstep+ coolmesh.array[0,0].x, -30 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) # top right
plt.scatter(350 * coolmesh.xstep+ coolmesh.array[0,0].x, -200 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) #bottome left
plt.scatter(650 * coolmesh.xstep+ coolmesh.array[0,0].x, -200 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) # bottom right
plt.show()

coolmesh.DumpFile("coolmesh2.msh")