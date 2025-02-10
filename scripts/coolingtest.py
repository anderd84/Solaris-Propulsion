import time
import matplotlib.pyplot as plt

from cooling import analysis, domain, domain_mmap, material
import general.design as DESIGN
from general.units import Q_, unitReg
from nozzle import plots, plug

def main():
    Re = Q_(3.2, unitReg.inch)
    exhaust = DESIGN.exhaustGas
    print(exhaust.stagTemp)

    cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
    Rt = outputData["radiusThroat"]
    Tt = outputData["thetaThroat"]
    Re = outputData["radiusLip"]

    # fig = plots.CreateNonDimPlot()
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

    # to run the saved one use this line:
    # coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh")
    coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("highmesh2")

    print(coolmesh.coolingLoops)

    mmapmesh = domain_mmap.DomainMMAP(coolmesh)

    fig = plots.CreateNonDimPlot()
    # print("speedrun to max diff 5%")
    tstart = time.time()
    analysis.AnalyzeMC(mmapmesh, 1, 1e-5, False, 5) # why wont you workkkkk 😭
    print(f"Time: {time.time() - tstart}")
    mmapmesh.NodePlot(fig, "material", [])
    plt.show()
    # analysis.AnalyzeMCSparse(coolmesh, 10, 1e-5, False)

if __name__ == "__main__":
    main()

