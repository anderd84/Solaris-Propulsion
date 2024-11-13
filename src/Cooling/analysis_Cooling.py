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

coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh.msh")

def getconductivity(coolmesh,i,j):
    match coolmesh.array[i,j].material:
        case DomainMaterial.COWL:
            return cooling_func.conduction_grcop(coolmesh.array[i,j].temperature.to(unitReg.degR))
        case DomainMaterial.COOLANT:
            return cooling_func.conduction_rp1(coolmesh.array[i,j].temperature.to(unitReg.degR))
        case DomainMaterial.PLUG:
            return cooling_func.conduction_grcop(coolmesh.array[i,j].temperature.to(unitReg.degR))                          



def getcorecool(coolmesh,i,j):
    #* If core & Coolant, Use the temperature below
    T_new = coolmesh.array[i+1,j].temperature
    return T_new

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

def coolant(coolmesh,i,j):
    """already on the boorder of coolant contour. first seperate alll top ones then go to the other ones


    Args:
        coolmesh (_type_): _description_
        i (_type_): _description_
        j (_type_): _description_

    Returns:
        _type_: _description_
    """

    if (coolmesh.array[i+1,j].material == DomainMaterial.COOLANT):
        T_new = coolmesh.array[i+1,j].temperature
    else:
    
        pass




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
            case DomainMaterial.COOLANT:
                conductance_conduction_left = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j-1].temperature.to(unitReg.degR),coolmesh.array[i,j-1].velocity.to(unitReg.foot/unitReg.second))
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
            case DomainMaterial.COOLANT:
                conductance_conduction_right = getconductivity(coolmesh,i,j) * (2 * np.pi * Q_(coolmesh.array[i,j].r, unitReg.inch).to(unitReg.foot) )* Q_(coolmesh.rstep, unitReg.inch).to(unitReg.foot) / (Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)/2)            
                convect = cooling_func.internal_flow_convection(coolmesh.array[i,j+1].temperature.to(unitReg.degR),coolmesh.array[i,j+1].velocity.to(unitReg.foot/unitReg.second))
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
                case DomainMaterial.COOLANT:
                    conduct_upper = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i, j].r  # Inner radius at [i, j] 
                    r_o = coolmesh.array[i-1,j].r - coolmesh.rstep/2  # Outer radius at [i-1,j] - halfstep
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_upper = (2 * np.pi * conduct_upper * l) / np.log(r_o / r_i)
                    convect = cooling_func.internal_flow_convection(coolmesh.array[i-1,j].temperature.to(unitReg.degR),coolmesh.array[i-1,j].velocity.to(unitReg.foot/unitReg.second))
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
                case DomainMaterial.COOLANT:
                    conduct_bottom = getconductivity(coolmesh,i,j)
                    r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                    r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                    l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                    conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                    convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),coolmesh.array[i+1,j].velocity.to(unitReg.foot/unitReg.second))
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
            case DomainMaterial.COOLANT:
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
            case DomainMaterial.COOLANT:
                conduct_bottom = getconductivity(coolmesh,i,j)
                r_i = coolmesh.array[i+1,j].r + coolmesh.rstep/2 # Inner radius at [i+1,j]  + half step
                r_o = coolmesh.array[i, j ].r   # Outer radius at [i, j]
                l = Q_(coolmesh.xstep, unitReg.inch).to(unitReg.foot)  # Axial length
                conductance_conduction_bottom = (2 * np.pi * conduct_bottom * l) / np.log(r_o / r_i)
                convect = cooling_func.internal_flow_convection(coolmesh.array[i+1,j].temperature.to(unitReg.degR),coolmesh.array[i+1,j].velocity.to(unitReg.foot/unitReg.second))
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

plt.ion()
plt.show()


for iterate in range(5):
    ic(iterate)
    fig.axes[0].clear()
    coolmesh.ShowStatePlot(fig)
    fig.canvas.draw()
    fig.canvas.flush_events()
    for i in range(coolmesh.vpoints):
        for j in range(coolmesh.hpoints):
            #*Finding all options for barrier
            if not(coolmesh.array[i,j].material == DomainMaterial.CHAMBER or coolmesh.array[i,j].material == DomainMaterial.FREE):

                if not(coolmesh.array[i,j].border):
                    if (coolmesh.array[i,j].material == DomainMaterial.COOLANT):
                        T_new = getcorecool(coolmesh,i,j)
                        coolmesh.array[i,j].temperature = Q_(T_new, unitReg.degR)
                        continue
                    else:
                        C_left, C_upper, C_bottom, C_right, T_left, T_upper, T_bottom, T_right = getcorecond(coolmesh,i,j)

                else:
                    if (coolmesh.array[i,j].material == DomainMaterial.COOLANT):
                        T_new = coolant(coolmesh,i,j)
                        coolmesh.array[i,j].temperature = Q_(T_new, unitReg.degR)
                        continue                                 
                    else:
                        C_left, T_left, C_right, T_right = horizontalcond(coolmesh,i,j)               
                        C_upper, T_upper, C_bottom, T_bottom = verticalcond(coolmesh,i,j)
                T_left, T_right, T_upper, T_bottom = Q_([T_left.magnitude, T_right.magnitude, T_upper.magnitude, T_bottom.magnitude], unitReg.degR)

                Num = (C_left * T_left + C_upper * T_upper + C_bottom * T_bottom + C_right*  T_right)
                Denom = (C_left + C_upper + C_bottom + C_right)
                coolmesh.array[i,j].temperature = (Num/Denom).to(unitReg.degR)



                
#TODO Add David's Gamma changing as a function of Temp





plt.show()