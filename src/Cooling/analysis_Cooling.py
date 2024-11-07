import matplotlib.pyplot as plt
from Cooling import domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg
from Cooling import cooling2d as cooling_func
from Cooling.material import DomainMaterial 


Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas
print(exhaust.stagTemp)

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

fig = plots.CreateNonDimPlot()
plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch), Q_(1.5, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
plots.PlotPlug(fig, chamberC, '-r')

coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh2.msh")


for i in range(coolmesh.vpoints):
    for j in range(coolmesh.hpoints):
        if coolmesh.array[i,j].border:
            if i==1:
                dasd = 2 #TODO placeholder fix 
            else:
                match coolmesh.array[i,j].material:
                    case COWL:
                        coolingfunc.
                       
                       
                        𝑇_𝑖=((𝑘 Δ𝑦/Δ𝑥 𝑇_𝐿+𝑘 Δ𝑦/Δ𝑥 𝑇_𝑅+" " 𝑘 Δ𝑥/Δ𝑦 𝑇_𝑇+𝑘 Δ𝑥/Δ𝑦 𝑇_𝐵 ))/((𝑘 Δ𝑦/Δ𝑥+𝑘 Δ𝑦/Δ𝑥+" " 𝑘 Δ𝑥/Δ𝑦+𝑘 Δ𝑥/Δ𝑦) )
                        
            

            if coolmesh.array[i,j].temperature:
                




#cooling_func.combustion_convection(mu, c_p, Pr, P_0, c_star, D_star, A_star, A, r_c, Ma, T_wg, T_0g, gamma)
#TODO Add David's Gamma changing as a function of Temp


coolmesh.ShowStatePlot(fig)



plt.show()