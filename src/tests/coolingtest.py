import matplotlib.pyplot as plt
from Cooling import analysis, analysisCoolingRef, domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg
from Cooling import cooling2d as cooling_func
from Cooling.material import DomainMaterial 
import numpy as np
from icecream import ic

def main():
    Re = Q_(3.2, unitReg.inch)
    exhaust = DESIGN.exhaustGas
    print(exhaust.stagTemp)

    cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
    Rt = outputData["radiusThroat"]
    Tt = outputData["thetaThroat"]
    Re = outputData["radiusLip"]

    # fig = plots.CreateNonDimPlot()
    plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
    cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch))
    chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.0203, unitReg.inch), Q_(1.5, unitReg.inch))
    # plots.PlotPlug(fig, plugC)
    # plots.PlotPlug(fig, cowlC)
    # plots.PlotPlug(fig, chamberC)
    # fig.axes[0].plot([p.x for p in cowlCoolL], [p.r for p in cowlCoolL], '-k', linewidth=1)
    # fig.axes[0].plot([p.x for p in cowlCoolU], [p.r for p in cowlCoolU], '-k', linewidth=1)
    # fig.axes[0].plot([p.x for p in plugCoolL], [p.r for p in plugCoolL], '-k', linewidth=1)
    # fig.axes[0].plot([p.x for p in plugCoolU], [p.r for p in plugCoolU], '-k', linewidth=1)

    # to run the saved one use this line:
    coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("save.msh")
    # coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh.msh")

    mmapmesh = domain.DomainMMAP(coolmesh)

    

    analysis.AnalyzeMC(mmapmesh, 0, plugC, cowlC)
    
    cool_mesh = mmapmesh.toDomain()
    cool_mesh.ShowStatePlot(0)
    cool_mesh.DumpFile("coolmesh.msh")
    plt.show()

if __name__ == "__main__":
    main()