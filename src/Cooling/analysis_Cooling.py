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





def leftnodeborder(coolmesh,i,j):
    #* Left Node
    if i==1:
        left_node_wo_T_L = Q_(0, unitReg.degR)
        deltaT_L = left_node_wo_T_L * coolmesh.array[i-1,j].temperature
    else:
        match coolmesh.array[i-1,j].material:
            case DomainMaterial.COWL:
                left_node_wo_T_L = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                deltaT_L = left_node_wo_T_L * coolmesh.array[i-1,j].temperature
            case DomainMaterial.CHAMBER:
                left_node_wo_T_L = cooling_func.conduction(coolmesh.array[i-1,j].temperature,coolmesh.array[i-1,j].velocity)
                deltaT_L = left_node_wo_T_L * coolmesh.array[i-1,j].temperature
    return left_node_wo_T_L, deltaT_L

def uppernodeborder(coolmesh,i,j):
    #* Upper Node
    if j==1:
        left_node_wo_T_U = Q_(0, unitReg.degR)
        deltaT_U = left_node_wo_T_U * coolmesh.array[i-1,j].temperature
    else:
        match coolmesh.array[i,j-1].material:
            case DomainMaterial.COWL:
                left_node_wo_T_U = cooling_func.conduction(coolmesh.array[i,j-1].temperature)
                deltaT_U = left_node_wo_T_U * coolmesh.array[i,j-1].temperature
            case DomainMaterial.CHAMBER:
                left_node_wo_T_U = cooling_func.conduction(coolmesh.array[i,j-1].temperature,coolmesh.array[i,j-1].velocity)
                deltaT_U = left_node_wo_T_U * coolmesh.array[i,j-1].temperature

    return left_node_wo_T_U, deltaT_U

def bottomnodeborder(coolmesh,i,j):
    #* Bottom Node
    if coolmesh.array[i,j] == coolmesh.array[-1,j]:
        left_node_wo_T_B = Q_(0, unitReg.degR)
        deltaT_B = left_node_wo_T_B * coolmesh.array[i+1,j].temperature
    else:
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                left_node_wo_T_B = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                deltaT_B = left_node_wo_T_B * coolmesh.array[i+1,j].temperature
            case DomainMaterial.CHAMBER:
                left_node_wo_T_B = cooling_func.conduction(coolmesh.array[i+1,j].temperature,coolmesh.array[i+1,j].velocity)ity)
                deltaT_B = left_node_wo_T_B * coolmesh.array[i+1,j].temperature
    return left_node_wo_T_B, deltaT_B

def rightnodeborder(coolmesh,i,j):
    #* Right Node
    if coolmesh.array[i,j] == coolmesh.array[i,-1]:
        left_node_wo_T_R = Q_(0, unitReg.degR)
        deltaT_R = left_node_wo_T_R * coolmesh.array[i,j+1].temperature
    else:
        match coolmesh.array[i,j+1].material:
            case DomainMaterial.COWL:
                left_node_wo_T_R = cooling_func.conduction(coolmesh.array[i,j+1].temperature)
                deltaT_R = left_node_wo_T_R * coolmesh.array[i,j+1].temperature
            case DomainMaterial.CHAMBER:
                left_node_wo_T_R = cooling_func.conduction(coolmesh.array[i,j+1].temperature,coolmesh.array[i,j+1].velocity)
                deltaT_R = left_node_wo_T_R * coolmesh.array[i,j+1].temperature
    return left_node_wo_T_R, deltaT_R


def calculate_nozzle_area(A_star, M, gamma):
    # Calculate area A using the rearranged equation
    A = A_star * ((gamma + 1) / 2) ** (- (gamma + 1) / (2 * (gamma - 1))) * (1 + (gamma - 1) / 2 * M**2) ** ((gamma + 1) / (2 * (gamma - 1))) * (1 / M)
    return A

for i in range(coolmesh.vpoints):
    for j in range(coolmesh.hpoints):
        #*Finding all options for barrier
        if coolmesh.array[i,j].border:
            left_node_wo_T_L, deltaT_L = leftnodeborder(coolmesh,i,j)
            left_node_wo_T_U, deltaT_U = uppernodeborder(coolmesh,i,j)
            left_node_wo_T_B, deltaT_B = bottomnodeborder(coolmesh,i,j)
            left_node_wo_T_R, deltaT_R = rightnodeborder(coolmesh,i,j)
            #* Finally calculate Temperature for Node
            Num = (deltaT_L + deltaT_U + deltaT_B + deltaT_R)
            Denom = (left_node_wo_T_L + left_node_wo_T_U + left_node_wo_T_B + left_node_wo_T_R)
            coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)


            
                
#TODO Add David's Gamma changing as a function of Temp


coolmesh.ShowStatePlot(fig)



plt.show()