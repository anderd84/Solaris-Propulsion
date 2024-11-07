import matplotlib.pyplot as plt
from Cooling import domain
from Nozzle import plots
import General.design as DESIGN
from Nozzle import plug
from General.units import Q_, unitReg
from Cooling import cooling2d as cooling_func
from Cooling.material import DomainMaterial 
import numpy as np


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
        left_node_wo_T = Q_(0, unitReg.degR)
        deltaT_L = left_node_wo_T * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i-1,j].material:
            case DomainMaterial.COWL:
                left_node_wo_T = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                deltaT_L = left_node_wo_T * coolmesh.array[i-1,j].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                convect = cooling_func.combustion_convection(coolmesh.array[i-1,j].temperature,coolmesh.array[i-1,j].velocity)
                
                left_node_wo_T = 1/(coolmesh.xstep/(conduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_L = left_node_wo_T * coolmesh.array[i-1,j].temperature
    return left_node_wo_T, deltaT_L

def uppernodeborder(coolmesh,i,j):
    #* Upper Node 
    if j==1:
        left_node_wo_T_U = Q_(0, unitReg.degR)
        deltaT_U = left_node_wo_T_U * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j + 1].material:
            case DomainMaterial.COWL:
                upper_node_wo_T = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Cylindrical conduction calculation using updated r_i and r_o
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                # Apply the cylindrical conduction formula
                upper_node_wo_T = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
    return upper_node_wo_T, deltaT_U

def bottomnodeborder(coolmesh,i,j):
    #* Bottom Node
    if coolmesh.array[i,j] == coolmesh.array[-1,j]:
        left_node_wo_T_B = Q_(0, unitReg.degR)
        deltaT_B = left_node_wo_T_B * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j - 1].material:
            case DomainMaterial.COWL:
                bottom_node_wo_T = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                # Cylindrical conduction calculation
                bottom_node_wo_T = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
    return bottom_node_wo_T, deltaT_B

def rightnodeborder(coolmesh,i,j):
    #* Right Node
    if coolmesh.array[i,j] == coolmesh.array[i,-1]:
        right_node_wo_T = Q_(0, unitReg.degR)
        deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature
    else:
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                right_node_wo_T = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                convect = cooling_func.combustion_convection(coolmesh.array[i+1,j].temperature,coolmesh.array[i+1,j].velocity)
                right_node_wo_T = 1/(coolmesh.xstep/(conduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature    
    return right_node_wo_T, deltaT_R


def calculate_nozzle_area(A_star, M, gamma):
    # Calculate area A using the rearranged equation
    A = A_star * ((gamma + 1) / 2) ** (- (gamma + 1) / (2 * (gamma - 1))) * (1 + (gamma - 1) / 2 * M**2) ** ((gamma + 1) / (2 * (gamma - 1))) * (1 / M)
    return A




for i in range(coolmesh.vpoints):
    for j in range(coolmesh.hpoints):
        #*Finding all options for barrier
        if coolmesh.array[i,j].border:
            match coolmesh.array[i,j].material:
                case DomainMaterial.COWL:
                    left_node_wo_T, deltaT_L = leftnodeborder(coolmesh,i,j)
                    upper_node_wo_T, deltaT_U = uppernodeborder(coolmesh,i,j)
                    bottom_node_wo_T, deltaT_B = bottomnodeborder(coolmesh,i,j)
                    right_node_wo_T, deltaT_R = rightnodeborder(coolmesh,i,j)
                    #* Finally calculate Temperature for Node
                    Num = (deltaT_L + deltaT_U + deltaT_B + deltaT_R)
                    Denom = (left_node_wo_T + upper_node_wo_T + bottom_node_wo_T + right_node_wo_T)
                    coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)
                case DomainMaterial.CHAMBER:
                    #DO nothing
                    


            
                
#TODO Add David's Gamma changing as a function of Temp


coolmesh.ShowStatePlot(fig)



plt.show()