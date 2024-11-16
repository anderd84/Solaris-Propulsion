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

check = coolmesh.array[60,529]

def getconductivity(coolmesh,i,j):
    match coolmesh.array[i,j].material:
        case DomainMaterial.COWL:
            return cooling_func.conduction_grcop(coolmesh.array[i,j].temperature.to(unitReg.degR))
        case DomainMaterial.COOLANT_WALL:
            return cooling_func.conduction_rp1(coolmesh.array[i,j].temperature.to(unitReg.degR))
        case DomainMaterial.COOLANT_BULK:
            return cooling_func.conduction_rp1(coolmesh.array[i,j].temperature.to(unitReg.degR))
        case DomainMaterial.PLUG:
            return cooling_func.conduction_grcop(coolmesh.array[i,j].temperature.to(unitReg.degR))                          


def getcorecond(coolmesh,i,j):
    #* Left Node
    if not(j==0):
        C_left = getconductivity(coolmesh,i,j-1) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) ) * Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
        T_left = coolmesh.array[i,j-1].temperature.to(unitReg.degR)
    else:
        C_left = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_left = Q_(0, unitReg.degR)
    #* Upper Node
    if not(i==0):
        conduct_upper = getconductivity(coolmesh,i-1,j)
        r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
        r_o = coolmesh.array[i - 1, j].r  # Outer radius at [i-1,j]
        l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
        cylindricalconduct = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
        C_upper = (cylindricalconduct) 
        T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)
    else:
        C_upper = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_upper = Q_(0, unitReg.degR)
    #* Bottom Node
    if not(coolmesh.array[i, j].r <=coolmesh.rstep):

        conduct_bottom = getconductivity(coolmesh,i+1,j)
        r_i = coolmesh.array[i + 1, j].r  # Inner radius at [i +1, j]
        r_o = coolmesh.array[i, j].r      # Outer radius at [i, j]
        l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)                # Axial length

        cylindricalconduct = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
        C_bottom = (cylindricalconduct)
        T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
    else:
        C_bottom = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_bottom = Q_(0, unitReg.degR)
    #* Right Node
    if not(j==coolmesh.hpoints-1):
        C_right = getconductivity(coolmesh,i,j + 1) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
        T_right = coolmesh.array[i,j + 1].temperature.to(unitReg.degR)
    else:
        C_right = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_right = Q_(0, unitReg.degR)

    return C_left, C_upper, C_bottom, C_right, T_left, T_upper, T_bottom, T_right

def coolant_wall_left(coolmesh,i,j):
    Area_exposed = (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) ) * Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)
    conduct_left = getconductivity(coolmesh,i,j-1) # gets the conductivity of the cell directly left of the previous node
    convect_left = cooling_func.internal_flow_convection((coolmesh.array[i,j].temperature).to(unitReg.degR),(coolmesh.array[i,j].pressure).to(unitReg.psi), (coolmesh.array[i,j].flowHeight).to(unitReg.inch))
    resistance_cond_left =  Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot )/ conduct_left * Area_exposed
    resistance_conv_left = 1 / convect_left * Area_exposed
    Q_left = ((coolmesh.array[i,j-1].temperature - coolmesh.array[i,j].temperature)/(resistance_cond_left+resistance_conv_left))
    Q_left = Q_left.to( unitReg.BTU / unitReg.hour)
    return Q_left

def coolant_wall_right(coolmesh,i,j):
    Area_exposed = (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) ) * Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)
    conduct_right = getconductivity(coolmesh,i,j+1) # gets the conductivity of the cell directly right of the previous node
    convect_right = cooling_func.internal_flow_convection((coolmesh.array[i,j].temperature).to(unitReg.degR),(coolmesh.array[i,j].pressure).to(unitReg.psi), (coolmesh.array[i,j].flowHeight).to(unitReg.inch))
    resistance_cond_right =  Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot )/ conduct_right * Area_exposed
    resistance_conv_right = 1 / convect_right * Area_exposed
    Q_right = (coolmesh.array[i,j+1].temperature - coolmesh.array[i,j].temperature)/(resistance_cond_right+resistance_conv_right)
    Q_right = Q_right.to( unitReg.BTU / unitReg.hour)
    return Q_right

def coolant_wall_below(coolmesh,i,j):
    Area_exposed = (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) ) * Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)
    conduct_below = getconductivity(coolmesh,i,j-1) # gets the conductivity of the cell directly below of the previous node
    r_i = coolmesh.array[i+1, j].r  # Inner radius at [i, j] 
    r_o = coolmesh.array[i+1,j].r + coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
    convect_below = cooling_func.internal_flow_convection((coolmesh.array[i,j].temperature).to(unitReg.degR),(coolmesh.array[i,j].pressure).to(unitReg.psi), (coolmesh.array[i,j].flowHeight).to(unitReg.inch))
    resistance_cond_below =  np.log(r_o / r_i) / (2 * np.pi * conduct_below * l)
    resistance_conv_below = 1 / convect_below * Area_exposed
    Q_below = (coolmesh.array[i+1,j].temperature - coolmesh.array[i,j].temperature)/(resistance_cond_below+resistance_conv_below)
    Q_below = Q_below.to( unitReg.BTU / unitReg.hour)
    return Q_below

def coolantwallheat(coolmesh,i,j, i_previous,j_previous):
    #*5 Cases overall
    if ((i_previous == i + 1) and (j_previous == j +1)):
    #*Case 1 (staircase)
        Q_left = coolant_wall_left(coolmesh,i_previous,j_previous)
        Q_below = coolant_wall_below(coolmesh,i_previous,j_previous)
        Qdotin = Q_left + Q_below #I need capital Q showings it's Heat rate in not Heat Flux #*WINSTON)
    elif ((i_previous == i) and (j_previous == j +1)):
    #*Case 2 & 3 (horizontal)

        if((coolmesh.array[i_previous,j_previous+1].material == DomainMaterial.PLUG) or(coolmesh.array[i_previous,j_previous+1].material == DomainMaterial.COWL)):
        #*Case 2 (2 qdots)
            Q_right = coolant_wall_right(coolmesh,i_previous,j_previous)
            Q_below = coolant_wall_below(coolmesh,i_previous,j_previous)
            Qdotin = Q_right + Q_below #I need capital Q showings it's Heat rate in not Heat Flux #*WINSTON)

        else:
        #*Case 3 (1 qdot1)
            Q_below = coolant_wall_below(coolmesh,i_previous,j_previous)
            Qdotin = Q_below #I need capital Q showings it's Heat rate in not Heat Flux #*WINSTON)

    elif ((i_previous == i+1) and (j_previous == j)):
    #*Case 4 & 5 (vertical)

        if((coolmesh.array[i_previous+1,j_previous].material == DomainMaterial.PLUG) or(coolmesh.array[i_previous+1,j_previous].material == DomainMaterial.COWL)):
        #*Case 4 (2 qdots)
            Q_left = coolant_wall_left(coolmesh,i_previous,j_previous)
            Q_below = coolant_wall_below(coolmesh,i_previous,j_previous)
            Qdotin = Q_left + Q_below #I need capital Q showings it's Heat rate in not Heat Flux #*WINSTON)

        else:
        #*Case 5 (1 qdot)
            Q_left = coolant_wall_left(coolmesh,i_previous,j_previous)
            Qdotin = Q_left #I need capital Q showings it's Heat rate in not Heat Flux #*WINSTON)
    else:   
        ic("YOu did something wrong for sure")
        Qdotin = 0
    return Qdotin


def coolant(coolmesh,i,j):


    if (coolmesh.array[i,j].material == DomainMaterial.COOLANT_WALL):
        i_previous = coolmesh.array[i,j].previousFlow[0]
        j_previous = coolmesh.array[i,j].previousFlow[1]
        Qdotin = coolantwallheat(coolmesh,i,j, i_previous,j_previous)
        T_new = cooling_func.heatcoolant(Qdotin, coolmesh.array[i_previous,j_previous].temperature.to(unitReg.degR), coolmesh.array[i_previous,j_previous].pressure)
    elif (coolmesh.array[i,j].material == DomainMaterial.COOLANT_BULK):
        i_wall = coolmesh.array[i,j].previousFlow[0]
        j_wall = coolmesh.array[i,j].previousFlow[1]
        T_new = coolmesh.array[i_wall,j_wall].temperature.to(unitReg.degR)
    return T_new

def horizontalcond(coolmesh,i,j):
    """ checks whether it is on the border of the mesh, then does vertical math to find conductances for 
        upper and lower nodes when the current cell is plug or cowl

    Args:
        coolmesh (_type_): _description_
        i (_type_): _description_
        j (_type_): _description_

    Returns:
        _type_: _description_
    """
    if not(j==0):
        match coolmesh.array[i,j-1].material:
            case DomainMaterial.COWL:
                C_left = getconductivity(coolmesh,i,j-1) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
                T_left = coolmesh.array[i,j-1].temperature.to(unitReg.degR)                
            case DomainMaterial.COOLANT_BULK:
                conductance_conduction_left = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j-1].temperature.to(unitReg.degR),(coolmesh.array[i,j-1].pressure).to(unitReg.psi), (coolmesh.array[i,j-1].flowHeight).to(unitReg.inch))
                conductance_convection_left = convect * (2 * np.pi * Q_(coolmesh.array[i, j].r, unitReg.inch).to(unitReg.foot) ) *Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)            
                C_left = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_left = coolmesh.array[i,j-1].temperature.to(unitReg.degR) 
            case DomainMaterial.COOLANT_WALL:
                conductance_conduction_left = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j-1].temperature.to(unitReg.degR),(coolmesh.array[i,j-1].pressure).to(unitReg.psi), (coolmesh.array[i,j-1].flowHeight).to(unitReg.inch))
                conductance_convection_left = convect * (2 * np.pi * Q_(coolmesh.array[i, j].r, unitReg.inch).to(unitReg.foot) ) *Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)            
                C_left = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_left = coolmesh.array[i,j-1].temperature.to(unitReg.degR)  
            case DomainMaterial.PLUG:
                C_left = getconductivity(coolmesh,i,j-1) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
                T_left = coolmesh.array[i,j-1].temperature.to(unitReg.degR)
            case DomainMaterial.CHAMBER:
                conductance_conduction_left = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.combustion_convection(coolmesh.array[i,j-1].temperature.to(unitReg.degR),coolmesh.array[i,j-1].velocity.to(unitReg.foot/unitReg.second))
                conductance_convection_left = convect * (2 * np.pi * Q_(coolmesh.array[i, j].r, unitReg.inch).to(unitReg.foot) ) * Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)            
                C_left = 1/(1/conductance_conduction_left + 1/conductance_convection_left)
                T_left = coolmesh.array[i,j-1].temperature.to(unitReg.degR)
            case _:
                C_left = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
                T_left = Q_(0, unitReg.degR)
    else:
        C_left = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_left = Q_(0, unitReg.degR)
    if not(j==coolmesh.hpoints-1):
        match coolmesh.array[i,j+1].material:
            case DomainMaterial.COWL:
                C_right = getconductivity(coolmesh,i,j+1) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
                T_right = coolmesh.array[i,j+1].temperature.to(unitReg.degR)               
            case DomainMaterial.COOLANT_BULK:
                conductance_conduction_right = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j+1].temperature.to(unitReg.degR),(coolmesh.array[i,j+1].pressure).to(unitReg.psi), (coolmesh.array[i,j+1].flowHeight).to(unitReg.inch))
                conductance_convection_right = convect *(2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)            
                C_right = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_right = coolmesh.array[i,j+1].temperature.to(unitReg.degR) 
            case DomainMaterial.COOLANT_WALL:
                conductance_conduction_right = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j+1].temperature.to(unitReg.degR),(coolmesh.array[i,j+1].pressure).to(unitReg.psi), (coolmesh.array[i,j+1].flowHeight).to(unitReg.inch))
                conductance_convection_right = convect *(2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)            
                C_right = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_right = coolmesh.array[i,j+1].temperature.to(unitReg.degR) 
            case DomainMaterial.PLUG:
                C_right = getconductivity(coolmesh,i,j+1) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
                T_right = coolmesh.array[i,j+1].temperature.to(unitReg.degR)
            case DomainMaterial.CHAMBER:
                conductance_conduction_right = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)      
                convect = cooling_func.combustion_convection(coolmesh.array[i,j+1].temperature.to(unitReg.degR),coolmesh.array[i,j+1].velocity.to(unitReg.foot/unitReg.second))
                conductance_convection_right = convect *(2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot)

                C_right = 1/(1/conductance_conduction_right + 1/conductance_convection_right)
                T_right = coolmesh.array[i,j+1].temperature.to(unitReg.degR)
            case _:
                C_right = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
                T_right = Q_(0, unitReg.degR)
    else:
        C_right =Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_right =Q_(0, unitReg.degR)
    return C_left, T_left, C_right, T_right


def verticalcond(coolmesh,i,j):
    """ checks whether it is on the border of the mesh, then does vertical math to find conductances for 
        upper and lower nodes when the current cell is plug or cowl

    Args:
        coolmesh (_type_): _description_
        i (_type_): _description_
        j (_type_): _description_

    Returns:
        _type_: _description_
    """    
    if not(coolmesh.array[i, j].r <=coolmesh.rstep):
        if not(i==0):
            match coolmesh.array[i-1,j].material:
                case DomainMaterial.COWL:
                    conduct_upper = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                    r_o = coolmesh.array[i-1,j].r  # Outer radius at [i-1,j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)

                    C_upper = conductance_conduction_upper
                    T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)              
                case DomainMaterial.COOLANT_BULK:
                    conduct_upper = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i, j].r  # Inner radius at [i, j] 
                    r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                    convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),(coolmesh.array[i-1,j].pressure).to(unitReg.psi), (coolmesh.array[i-1,j].flowHeight).to(unitReg.inch))
                    conductance_convection_upper = convect * (2 * np.pi * Q_( coolmesh.array[i-1,j].r - coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                    C_upper = 1/(1/conductance_conduction_upper + 1/conductance_convection_upper)
                    T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)
                case DomainMaterial.COOLANT_WALL:
                    conduct_upper = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i, j].r  # Inner radius at [i, j] 
                    r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                    convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),(coolmesh.array[i-1,j].pressure).to(unitReg.psi), (coolmesh.array[i-1,j].flowHeight).to(unitReg.inch))
                    conductance_convection_upper = convect * (2 * np.pi * Q_( coolmesh.array[i-1,j].r - coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                    C_upper = 1/(1/conductance_conduction_upper + 1/conductance_convection_upper)
                    T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)                
                case DomainMaterial.PLUG:
                    conduct_upper = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i, j].r  # Inner radius at [i, j]
                    r_o = coolmesh.array[i-1,j].r  # Outer radius at [i-1,j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                    C_upper = conductance_conduction_upper
                    T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)
                case DomainMaterial.CHAMBER:
                    conduct_upper = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i, j].r  # Inner radius at [i, j] 
                    r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                    convect = cooling_func.combustion_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),coolmesh.array[i-1,j].velocity.to(unitReg.foot/unitReg.second))
                    conductance_convection_upper = convect * (2 * np.pi * Q_( coolmesh.array[i-1,j].r - coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                    C_upper = 1/(1/conductance_conduction_upper + 1/conductance_convection_upper)
                    T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR) 
                case _:
                    C_upper = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
                    T_upper = Q_(0, unitReg.degR)
            match coolmesh.array[i+1,j].material:
                case DomainMaterial.COWL:
                    conduct_bottom = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i+1,j].r  # Inner radius at [i+1,j]
                    r_o = coolmesh.array[i, j].r  # Outer radius at [i, j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                    C_bottom = conductance_conduction_bottom
                    T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)                
                case DomainMaterial.COOLANT_BULK:
                    conduct_bottom = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                    r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                    convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),(coolmesh.array[i+1,j].pressure).to(unitReg.psi), (coolmesh.array[i+1,j].flowHeight).to(unitReg.inch))
                    conductance_convection_bottom = convect * (2 * np.pi * Q_( coolmesh.array[i+1,j].r + coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                    C_bottom = 1/(1/conductance_conduction_bottom + 1/conductance_convection_bottom)
                    T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
                case DomainMaterial.COOLANT_WALL:
                    conduct_bottom = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                    r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                    convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),(coolmesh.array[i+1,j].pressure).to(unitReg.psi), (coolmesh.array[i+1,j].flowHeight).to(unitReg.inch))
                    conductance_convection_bottom = convect * (2 * np.pi * Q_( coolmesh.array[i+1,j].r + coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                    C_bottom = 1/(1/conductance_conduction_bottom + 1/conductance_convection_bottom)
                    T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
                case DomainMaterial.PLUG:
                    conduct_bottom = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i+1,j].r  # Inner radius at [i, j-1]
                    r_o = coolmesh.array[i, j].r  # Outer radius at [i, j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_bottom =  (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                    C_bottom = conductance_conduction_bottom
                    T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
                case DomainMaterial.CHAMBER:
                    conduct_bottom = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                    r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                    convect = cooling_func.combustion_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),coolmesh.array[i+1,j].velocity.to(unitReg.foot/unitReg.second))
                    conductance_convection_bottom = convect * (2 * np.pi * Q_( coolmesh.array[i+1,j].r + coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
                    C_bottom = 1/(1/conductance_conduction_bottom + 1/conductance_convection_bottom)
                    T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
                case _:
                    C_bottom = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
                    T_bottom = Q_(0, unitReg.degR)
        else:
            C_upper =Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
            T_upper =Q_(0, unitReg.degR)
            C_bottom = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
            T_bottom = Q_(0, unitReg.degR)  
    elif (coolmesh.array[i, j].r == coolmesh.rstep):
        match coolmesh.array[i-1,j].material:
            case DomainMaterial.COWL:
                conduct_upper = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10
                r_o = coolmesh.array[i-1,j].r  # Outer radius at [i-1,j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)

                C_upper = conductance_conduction_upper
                T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)              
            case DomainMaterial.COOLANT_BULK:
                conduct_upper = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10
                r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),(coolmesh.array[i-1,j].pressure).to(unitReg.psi), (coolmesh.array[i-1,j].flowHeight).to(unitReg.inch))
                conductance_convection_upper = convect * (2 * np.pi * Q_( coolmesh.array[i-1,j].r - coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                C_upper = 1/(1/conductance_conduction_upper + 1/conductance_convection_upper)
                T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)
            case DomainMaterial.COOLANT_WALL:
                conduct_upper = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10
                r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),coolmesh.array[i-1,j].velocity.to(unitReg.foot/unitReg.second))
                conductance_convection_upper = convect * (2 * np.pi * Q_( coolmesh.array[i-1,j].r - coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                C_upper = 1/(1/conductance_conduction_upper + 1/conductance_convection_upper)
                T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)
            case DomainMaterial.PLUG:
                conduct_upper = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10
                r_o = coolmesh.array[i-1,j].r  # Outer radius at [i-1,j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                C_upper = conductance_conduction_upper
                T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR)
            case DomainMaterial.CHAMBER:
                conduct_upper = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10 
                r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                convect = cooling_func.combustion_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),coolmesh.array[i-1,j].velocity.to(unitReg.foot/unitReg.second))
                conductance_convection_upper = convect * (2 * np.pi * Q_( coolmesh.array[i-1,j].r - coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                C_upper = 1/(1/conductance_conduction_upper + 1/conductance_convection_upper)
                T_upper = coolmesh.array[i-1,j].temperature.to(unitReg.degR) 
            case _:
                C_upper = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
                T_upper = Q_(0, unitReg.degR)
        match coolmesh.array[i+1,j].material:
            case DomainMaterial.COWL:
                conduct_bottom = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10
                r_o = coolmesh.array[i, j].r  # Outer radius at [i, j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                C_bottom = conductance_conduction_bottom
                T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)                
            case DomainMaterial.COOLANT_BULK:
                conduct_bottom = getconductivity(coolmesh,i,j)
                r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),(coolmesh.array[i+1,j].pressure).to(unitReg.psi), (coolmesh.array[i+1,j].flowHeight).to(unitReg.inch))
                conductance_convection_bottom = convect * (2 * np.pi * Q_( coolmesh.array[i+1,j].r + coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                C_bottom = 1/(1/conductance_conduction_bottom + 1/conductance_convection_bottom)
                T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
            case DomainMaterial.COOLANT_WALL:
                conduct_bottom = getconductivity(coolmesh,i,j)
                r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),(coolmesh.array[i+1,j].pressure).to(unitReg.psi), (coolmesh.array[i+1,j].flowHeight).to(unitReg.inch))
                conductance_convection_bottom = convect * (2 * np.pi * Q_( coolmesh.array[i+1,j].r + coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)            
                C_bottom = 1/(1/conductance_conduction_bottom + 1/conductance_convection_bottom)
                T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
            case DomainMaterial.PLUG:
                conduct_bottom = getconductivity(coolmesh,i,j)
                r_i = coolmesh.rstep/10
                r_o = coolmesh.array[i, j].r  # Outer radius at [i, j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_bottom =  (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                C_bottom = conductance_conduction_bottom
                T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
            case DomainMaterial.CHAMBER:
                conduct_bottom = getconductivity(coolmesh,i,j)
                r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                convect = cooling_func.combustion_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),coolmesh.array[i+1,j].velocity.to(unitReg.foot/unitReg.second))
                conductance_convection_bottom = convect * (2 * np.pi * Q_( coolmesh.array[i+1,j].r + coolmesh.rstep/2, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)
                C_bottom = 1/(1/conductance_conduction_bottom + 1/conductance_convection_bottom)

                T_bottom = coolmesh.array[i+1,j].temperature.to(unitReg.degR)
            case _:
                C_bottom = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
                T_bottom = Q_(0, unitReg.degR)
    else:
        C_upper =Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_upper =Q_(0, unitReg.degR)
        C_bottom = Q_(0, unitReg.BTU  / unitReg.hour / unitReg.degR)
        T_bottom = Q_(0, unitReg.degR)   
    
    return C_upper, T_upper, C_bottom, T_bottom

#plt.ion()
#plt.show()

ic(coolmesh.array[59,350].border)
ic(coolmesh.array[59,350].material)
ic(coolmesh.array[60,350].material)
ic(coolmesh.array[58,350].material)
ic(coolmesh.array[59,351].material)
ic(coolmesh.array[59,349].material)

#
#

total_change = 0
for iterate in range(5):
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
            if not(coolmesh.array[i,j].material == DomainMaterial.CHAMBER or coolmesh.array[i,j].material == DomainMaterial.FREE):  # Select only walls and coolant
                
                if (coolmesh.array[i,j].material == DomainMaterial.COOLANT_WALL or coolmesh.array[i,j].material == DomainMaterial.COOLANT_BULK):    # Select coolant
                        T_new = coolant(coolmesh,i,j)
                        coolmesh.array[i,j].temperature = Q_(T_new.magnitude, unitReg.degR)
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
            



                
#TODO Add David's Gamma changing as a function of Temp
coolmesh.ShowStatePlot(fig)
plt.scatter(350 * coolmesh.xstep+ coolmesh.array[0,0].x, -30 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) #top left
plt.scatter(650 * coolmesh.xstep+ coolmesh.array[0,0].x, -30 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) # top right
plt.scatter(350 * coolmesh.xstep+ coolmesh.array[0,0].x, -200 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) #bottome left
plt.scatter(650 * coolmesh.xstep+ coolmesh.array[0,0].x, -200 *coolmesh.rstep + coolmesh.array[0,0].r, marker='x', s=100) # bottom right
plt.show()