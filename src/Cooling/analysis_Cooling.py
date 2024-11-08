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
plugC, straightLength = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(5, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch))
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(5, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, Q_(0.025, unitReg.inch), Q_(1.5, unitReg.inch))
plots.PlotPlug(fig, plugC)
plots.PlotPlug(fig, cowlC)
plots.PlotPlug(fig, chamberC, '-r')

coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh2.msh")

def getconductivity(coolmesh,i,j):
    match coolmesh.array[i,j].material:
        case DomainMaterial.COWL:
            cooling_func.conduction_grcop(coolmesh.array[i,j].temperature)
        case DomainMaterial.COOLANT:
            cooling_func.conduction_rp1(coolmesh.array[i,j].temperature)
        case DomainMaterial.PLUG:
            cooling_func.conduction_grcop(coolmesh.array[i,j].temperature)                          

def getcore(coolmesh,i,j):
    #* Left Node
    if not(i==0) :
        C_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / coolmesh.xstep
        T_left = coolmesh.array[i-1,j].temperature
    else:
        C_left =0
        T_left =0
    #* Upper Node
    if not(j==0):
        conduct_upper = getconductivity(coolmesh,i,j+1)
        r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
        r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
        l = coolmesh.xstep  # Axial length
        cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct_upper * l)
        C_upper = (1/cylindricalconduct) 
        T_upper = coolmesh.array[i,j-1].temperature
    else:
        C_upper = 0
        T_upper = 0
    #* Bottom Node
    if not(j==coolmesh.hpoints-1):

        conduct_bottom = getconductivity(coolmesh,i,j-1)
        r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
        r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
        l = coolmesh.xstep                # Axial length

        cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct_bottom * l)
        C_bottom = (1/cylindricalconduct)
        T_bottom = coolmesh.array[i,j-1].temperature
    else:
        C_bottom = 0
        T_bottom = 0
    #* Right Node
    if not(i==coolmesh.vpoints-1):
        C_right = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / coolmesh.xstep
        T_right = coolmesh.array[i-1,j].temperature
    else:
        C_right = 0
        T_right = 0

   
    return C_left, C_upper, C_bottom, C_right, T_left, T_upper, T_bottom, T_right

def horizontalcool(coolmesh,i,j):
    if not(i==0) :
        match coolmesh.array[i-1,j].material:
            case DomainMaterial.COWL:
                conductance_conduction_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j].temperature,coolmesh.array[i,j].velocity)
                conductance_convection_left = convect * coolmesh.rstep            
                C_left = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_left = coolmesh.array[i-1,j].temperature                
            case DomainMaterial.COOLANT:
                C_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / coolmesh.xstep
                T_left = coolmesh.array[i-1,j].temperature

    else:
        C_left =0
        T_left =0
    if not(i==coolmesh.vpoints-1) :
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                conductance_conduction_right = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j].temperature,coolmesh.array[i,j].velocity)
                conductance_convection_right = convect * coolmesh.rstep            
                C_right = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_right = coolmesh.array[i+1,j].temperature                
            case DomainMaterial.COOLANT:
                C_right = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / coolmesh.xstep
                T_right = coolmesh.array[i+1,j].temperature
 
    else:
        C_right =0
        T_right =0    
    
    return C_left, T_left, C_right, T_right

def horizontalcond(coolmesh,i,j):
    if not(i==0) :
        match coolmesh.array[i-1,j].material:
            case DomainMaterial.COWL:
                C_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / coolmesh.xstep
                T_left = coolmesh.array[i-1,j].temperaturee                
            case DomainMaterial.COOLANT:
                conductance_conduction_left = getconductivity(coolmesh,i,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature,coolmesh.array[i-1,j].velocity)
                conductance_convection_left = convect * coolmesh.rstep            
                C_left = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_left = coolmesh.array[i-1,j].temperature  
            case DomainMaterial.PLUG:
                C_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / coolmesh.xstep
                T_left = coolmesh.array[i-1,j].temperature
            case DomainMaterial.CHAMBER:
                conductance_conduction_left = getconductivity(coolmesh,i,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i-1,j].temperature,coolmesh.array[i-1,j].velocity)
                conductance_convection_left = convect * coolmesh.rstep            
                C_left = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_left = coolmesh.array[i-1,j].temperature 
    else:
        C_left =0
        T_left =0
    if not(i==coolmesh.vpoints-1) :
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                C_right = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / coolmesh.xstep
                T_right = coolmesh.array[i+1,j].temperature               
            case DomainMaterial.COOLANT:
                conductance_conduction_right = getconductivity(coolmesh,i,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature,coolmesh.array[i+1,j].velocity)
                conductance_convection_right = convect * coolmesh.rstep            
                C_right = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_right = coolmesh.array[i+1,j].temperature 
            case DomainMaterial.PLUG:
                C_right = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / coolmesh.xstep
                T_right = coolmesh.array[i+1,j].temperature
            case DomainMaterial.CHAMBER:
                conductance_conduction_right = getconductivity(coolmesh,i,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i+1,j].temperature,coolmesh.array[i+1,j].velocity)
                conductance_convection_right = convect * coolmesh.rstep            
                C_right = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_right = coolmesh.array[i+1,j].temperature 
    else:
        C_right =0
        T_right =0    
    
    return C_left, T_left, C_right, T_right

def verticalcool(coolmesh,i,j):
    if not(j==0) :
        match coolmesh.array[i,j+1].material:
            case DomainMaterial.COWL:
                conductance_conduction_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i,j].temperature,coolmesh.array[i,j].velocity)
                conductance_convection_left = convect * coolmesh.rstep            
                C_upper = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_upper = coolmesh.array[i-1,j].temperature                
            case DomainMaterial.COOLANT:
                C_upper = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / coolmesh.xstep
                T_upper = coolmesh.array[i-1,j].temperature
            case DomainMaterial.PLUG:
                conductance_conduction_left = getconductivity(coolmesh,i-1,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i,j].temperature,coolmesh.array[i,j].velocity)
                conductance_convection_left = convect * coolmesh.rstep            
                C_upper = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_upper = coolmesh.array[i-1,j].temperature 
    else:
        C_upper =0
        T_upper =0
    if not(j==coolmesh.hpoints-1) :
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                conductance_conduction_right = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i,j].temperature,coolmesh.array[i,j].velocity)
                conductance_convection_right = convect * coolmesh.rstep            
                C_bottom = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_bottom = coolmesh.array[i+1,j].temperature                
            case DomainMaterial.COOLANT:
                C_bottom = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / coolmesh.xstep
                T_bottom = coolmesh.array[i+1,j].temperature
            case DomainMaterial.PLUG:
                conductance_conduction_right = getconductivity(coolmesh,i+1,j) * coolmesh.rstep / (coolmesh.xstep/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i,j].temperature,coolmesh.array[i,j].velocity)
                conductance_convection_right = convect * coolmesh.rstep            
                C_bottom = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_bottom = coolmesh.array[i+1,j].temperature 
    else:
        C_bottom =0
        T_bottom =0  
    

    return C_upper, T_upper, C_bottom, T_bottom






    #* Upper Node 
    if j==0:
        left_node_wo_T_U = Q_(0, unitReg.degR)
        deltaT_U = left_node_wo_T_U * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j + 1].material:
            case DomainMaterial.COWL:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                convect = cooling_func.combustion_convection(coolmesh.array[i,j+1].temperature,coolmesh.array[i,j+1].velocity)           
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.COOLANT:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j+1].temperature,coolmesh.array[i,j+1].velocity)           
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.PLUG:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
























def verticalcond(coolmesh,i,j):
    pass #TODO still doing the thing not actuall pass
    return C_bottom, T_bottom, C_right, T_right

def getleftright(coolmesh,i,j):    
    if (coolmesh.array[i,j].material == DomainMaterial.COOLANT):
        C_left, T_left, C_right, T_right = horizontalcool(coolmesh,i,j)
    else:
        C_left, T_left, C_right, T_right = horizontalcond(coolmesh,i,j)    
    return C_left, T_left, C_right, T_right
            

def getupperbottom(coolmesh,i,j):           
    if (coolmesh.array[i,j].material == DomainMaterial.COOLANT):
        C_upper, T_upper, C_bottom, T_bottom = verticalcool(coolmesh,i,j)
    else:
        C_upper, T_upper, C_bottom, T_bottom = verticalcond(coolmesh,i,j)    
    return C_upper, T_upper, C_bottom, T_bottom


for i in range(coolmesh.vpoints):
    for j in range(coolmesh.hpoints):
        #*Finding all options for barrier
        if not(coolmesh.array[i,j].material == DomainMaterial.CHAMBER):
            if not(coolmesh.array[i,j].border):
                C_left, C_upper, C_bottom, C_right, T_left, T_upper, T_bottom, T_right = getcore(coolmesh,i,j)

            else:
                C_left, T_left, C_right, T_right = getleftright(coolmesh,i,j)                
                C_upper, T_upper, C_bottom, T_bottom = getupperbottom(coolmesh,i,j)                        

            Num = (C_left * T_left + C_upper * T_upper + C_bottom * T_bottom + C_right*  T_right)
            Denom = (C_left + C_upper + C_bottom + C_right)
            coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)
        



















#*Below this is the old version that I made (never been tested)

def con_leftnodeborder(coolmesh,i,j):
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
            case DomainMaterial.COOLANT:
                conduct = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature,coolmesh.array[i-1,j].velocity)
                left_node_wo_T = 1/(coolmesh.xstep/(conduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_L = left_node_wo_T * coolmesh.array[i-1,j].temperature    
            case DomainMaterial.PLUG:
                left_node_wo_T = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                deltaT_L = left_node_wo_T * coolmesh.array[i-1,j].temperature                        
    return left_node_wo_T, deltaT_L

def con_uppernodeborder(coolmesh,i,j):
    #* Upper Node 
    if j==1:
        left_node_wo_T_U = Q_(0, unitReg.degR)
        deltaT_U = left_node_wo_T_U * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j + 1].material:
            case DomainMaterial.COWL:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                convect = cooling_func.combustion_convection(coolmesh.array[i,j+1].temperature,coolmesh.array[i,j+1].velocity)           
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.COOLANT:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j+1].temperature,coolmesh.array[i,j+1].velocity)           
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
            case DomainMaterial.PLUG:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
    return upper_node_wo_T, deltaT_U

def con_bottomnodeborder(coolmesh,i,j):
    #* Bottom Node
    if coolmesh.array[i,j] == coolmesh.array[-1,j]:
        left_node_wo_T_B = Q_(0, unitReg.degR)
        deltaT_B = left_node_wo_T_B * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j - 1].material:
            case DomainMaterial.COWL:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                bottom_node_wo_T = 1 / (coolmesh.xstep / (cylindricalconduct * coolmesh.rstep))
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
            case DomainMaterial.CHAMBER:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                convect = cooling_func.combustion_convection(coolmesh.array[i, j - 1].temperature, coolmesh.array[i, j - 1].velocity)
                bottom_node_wo_T = 1 / (coolmesh.xstep / (cylindricalconduct * coolmesh.rstep) + 1 / (convect * coolmesh.rstep))
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
            case DomainMaterial.COOLANT:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i, j - 1].temperature, coolmesh.array[i, j - 1].velocity)
                bottom_node_wo_T = 1 / (coolmesh.xstep / (cylindricalconduct * coolmesh.rstep) + 1 / (convect * coolmesh.rstep))
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
            case DomainMaterial.PLUG:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                bottom_node_wo_T = 1 / (coolmesh.xstep / (cylindricalconduct * coolmesh.rstep))
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
    return bottom_node_wo_T, deltaT_B

def con_rightnodeborder(coolmesh,i,j):
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
            case DomainMaterial.COOLANT:
                conduct = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature,coolmesh.array[i+1,j].velocity)
                right_node_wo_T = 1/(coolmesh.xstep/(conduct * coolmesh.rstep) + 1/(convect * coolmesh.rstep))
                deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature
            case DomainMaterial.PLUG:
                right_node_wo_T = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature         
    return right_node_wo_T, deltaT_R



def con_leftnodecore(coolmesh,i,j):
    #* Left Node
    if i==1:
        left_node_wo_T = Q_(0, unitReg.degR)
        deltaT_L = left_node_wo_T * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i-1,j].material:
            case DomainMaterial.COWL:
                left_node_wo_T = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                deltaT_L = left_node_wo_T * coolmesh.array[i-1,j].temperature
            case DomainMaterial.PLUG:
                left_node_wo_T = cooling_func.conduction(coolmesh.array[i-1,j].temperature)
                deltaT_L = left_node_wo_T * coolmesh.array[i-1,j].temperature                        
    return left_node_wo_T, deltaT_L

def con_uppernodecore(coolmesh,i,j):
    #* Upper Node 
    if j==1:
        left_node_wo_T_U = Q_(0, unitReg.degR)
        deltaT_U = left_node_wo_T_U * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j + 1].material:
            case DomainMaterial.COWL:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature

            case DomainMaterial.PLUG:
                conduct = cooling_func.conduction(coolmesh.array[i, j + 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                r_o = coolmesh.array[i, j + 1].r  # Outer radius at [i, j+1]
                l = coolmesh.xstep  # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                upper_node_wo_T = 1/(coolmesh.xstep/(cylindricalconduct * coolmesh.rstep))
                deltaT_U = upper_node_wo_T * coolmesh.array[i, j + 1].temperature
    return upper_node_wo_T, deltaT_U

def con_bottomnodecore(coolmesh,i,j):
    #* Bottom Node
    if coolmesh.array[i,j] == coolmesh.array[-1,j]:
        left_node_wo_T_B = Q_(0, unitReg.degR)
        deltaT_B = left_node_wo_T_B * coolmesh.array[i,j].temperature
    else:
        match coolmesh.array[i, j - 1].material:
            case DomainMaterial.COWL:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                bottom_node_wo_T = 1 / (coolmesh.xstep / (cylindricalconduct * coolmesh.rstep))
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
            case DomainMaterial.PLUG:
                conduct = cooling_func.conduction(coolmesh.array[i, j - 1].temperature)
                
                # Define inner and outer radii for cylindrical conduction
                r_i = coolmesh.array[i, j - 1].r  # Inner radius at [i, j-1]
                r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
                l = coolmesh.xstep                # Axial length

                cylindricalconduct = np.log(r_o / r_i) / (2 * np.pi * conduct * l)
                bottom_node_wo_T = 1 / (coolmesh.xstep / (cylindricalconduct * coolmesh.rstep))
                deltaT_B = bottom_node_wo_T * coolmesh.array[i, j - 1].temperature
    return bottom_node_wo_T, deltaT_B

def con_rightnodecore(coolmesh,i,j):
    #* Right Node
    if coolmesh.array[i,j] == coolmesh.array[i,-1]:
        right_node_wo_T = Q_(0, unitReg.degR)
        deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature
    else:
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                right_node_wo_T = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature
            case DomainMaterial.PLUG:
                right_node_wo_T = cooling_func.conduction(coolmesh.array[i+1,j].temperature)
                deltaT_R = right_node_wo_T * coolmesh.array[i+1,j].temperature         
    return right_node_wo_T, deltaT_R


for i in range(coolmesh.vpoints):
    for j in range(coolmesh.hpoints):
        #*Finding all options for barrier
        if coolmesh.array[i,j].border:
            match coolmesh.array[i,j].material:
                case DomainMaterial.COWL:
                    left_node_wo_T, deltaT_L = con_leftnodeborder(coolmesh,i,j)
                    upper_node_wo_T, deltaT_U = con_uppernodeborder(coolmesh,i,j)
                    bottom_node_wo_T, deltaT_B = con_bottomnodeborder(coolmesh,i,j)
                    right_node_wo_T, deltaT_R = con_rightnodeborder(coolmesh,i,j)
                    #* Finally calculate Temperature for Node
                    Num = (deltaT_L + deltaT_U + deltaT_B + deltaT_R)
                    Denom = (left_node_wo_T + upper_node_wo_T + bottom_node_wo_T + right_node_wo_T)
                    coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)

                case DomainMaterial.CHAMBER:
                    pass

                case DomainMaterial.PLUG:
                    left_node_wo_T, deltaT_L = con_leftnodeborder(coolmesh,i,j)
                    upper_node_wo_T, deltaT_U = con_uppernodeborder(coolmesh,i,j)
                    bottom_node_wo_T, deltaT_B = con_bottomnodeborder(coolmesh,i,j)
                    right_node_wo_T, deltaT_R = con_rightnodeborder(coolmesh,i,j)
                    #* Finally calculate Temperature for Node
                    Num = (deltaT_L + deltaT_U + deltaT_B + deltaT_R)
                    Denom = (left_node_wo_T + upper_node_wo_T + bottom_node_wo_T + right_node_wo_T)
                    coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)

        else:

            match coolmesh.array[i,j].material:
                case DomainMaterial.COWL:
                    left_node_wo_T, deltaT_L = con_leftnodecore(coolmesh,i,j)
                    upper_node_wo_T, deltaT_U = con_uppernodecore(coolmesh,i,j)
                    bottom_node_wo_T, deltaT_B = con_bottomnodecore(coolmesh,i,j)
                    right_node_wo_T, deltaT_R = con_rightnodecore(coolmesh,i,j)
                    #* Finally calculate Temperature for Node
                    Num = (deltaT_L + deltaT_U + deltaT_B + deltaT_R)
                    Denom = (left_node_wo_T + upper_node_wo_T + bottom_node_wo_T + right_node_wo_T)
                    coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)

                case DomainMaterial.CHAMBER:
                    a = 0 #doing nothing
                    
                case DomainMaterial.PLUG:
                    left_node_wo_T, deltaT_L = con_leftnodeborder(coolmesh,i,j)
                    upper_node_wo_T, deltaT_U = con_uppernodeborder(coolmesh,i,j)
                    bottom_node_wo_T, deltaT_B = con_bottomnodeborder(coolmesh,i,j)
                    right_node_wo_T, deltaT_R = con_rightnodeborder(coolmesh,i,j)
                    #* Finally calculate Temperature for Node
                    Num = (deltaT_L + deltaT_U + deltaT_B + deltaT_R)
                    Denom = (left_node_wo_T + upper_node_wo_T + bottom_node_wo_T + right_node_wo_T)
                    coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)




            
                
#TODO Add David's Gamma changing as a function of Temp


coolmesh.ShowStatePlot(fig)



plt.show()