from matplotlib import pyplot as plt, patches
from Doublet_Functions import spike_contour
from InjectorCad import injector_cad_write
from Drill import drill_approximation
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from pint import UnitRegistry
import numpy as np
import os
import pyromat as pm
import subprocess

"""
    Shit that needs to get done still in this code:
    Utilize the Space prop notes to find time of vaporization and make sure our chamber length design fits for that. For instance the larger the diameter
    the larger the time necessary to vaporize and the larger the L* we need
    Add the film cooling angles onto the plot
    have A seperate figure for a 2D sketch of the hgoles to make visualization much much easier
    Error handle the numerical solve for if the angle is bigger than 90 or less than the minimum angle possible (USE TRY CATCH)
    Make all code blocks into functions - EASIER TO READ. So david no kill me
"""

pm.config['unit_pressure'] = 'psi'
pm.config['unit_mass'] = 'lbm'
pm.config['unit_matter'] = 'lbm'
pm.config['unit_length'] = 'ft'
pm.config['unit_volume'] = 'ft3'
ureg = UnitRegistry()
ureg.default_system = 'US'
Q_ = ureg.Quantity

N2 = pm.get('mp.N2')
LOX = pm.get('mp.O2')
CO2 = pm.get('mp.CO2')

class PROP:
    def __init__(self, gamma, mdot, rho):
        """This is teh __init__ fn that will do a thing
        Args:
            gamma (float): Angle off axial direction with positive away from centerbody
            mdot (float): mdot of prop
            rho (float): density of prop
        """        
        self.gamma = Q_(gamma, ureg.degrees)
        self.mdot = Q_(mdot, ureg.pound / ureg.second)
        self.rho = Q_(rho, ureg.pound / ureg.foot**3)


    def Velocity(self,Cd , Pressure_Diff) -> float:
        """Velocity Function using Dr. Whites Eq 16.5
        Args:
            Cd (float): Coefficient of Discharge
            Pressure_Diff (float): Pressure difference through injector (Estimated metric) in psi
        Returns:
            float: Velocity Out of the injector parallel to flow
        """
        Velocity_init = (Cd * np.sqrt(2 * Pressure_Diff/ self.rho )).to(ureg.feet / ureg.second)
        self.Velocity_init = Velocity_init
        return Velocity_init
    
    
    def Area(self, Cd, Pressure_Diff) -> float:
        """Calculates area from Dr. Whites Eq 16.2
        Args:
            Cd (float): Coefficient of Discharge
            Pressure_Diff (float): Pressure difference through injector (Estimated metric) in psi
        Returns:
            float: Total Orifice Area needed for Propellant
        """        
        Area_init = (self.mdot / (Cd * np.sqrt(2*self.rho* Pressure_Diff ))).to(ureg.inch**2)
        self.Area_init = Area_init
        return Area_init
    
    
    def Number(self, Hole_Diameter, Cd, Pressure_Diff) -> float:
        """_summary_
        Args:
            Hole_Diameter (float): Diameter necessary for hole (MUST MATCH A DRILL SIZE)
            Cd (float): Coefficient of Discharge
            Pressure_Diff (float): Pressure difference through injector (Estimated metric) in psi
        Returns:
            float: Number of holes needed rounded up to be conservative
        """        
        Hole_Area = (np.pi * Hole_Diameter**2 /4).to(ureg.inch**2)
        Tot_Area = self.Area(Cd,Pressure_Diff)
        Hole_Number = np.ceil(Tot_Area/Hole_Area)
        self.Hole_Number = Hole_Number
        return Hole_Number
    def Actual(self,Hole_Diameter, Number) -> float:
        """_summary_

        Args:
            Hole_Diameter (float): Diameter of the holes, WIll match a drill size
            Number (float): Calculated number of holes that matched the other impingement for core impingement points

        Returns:
            float: Returns actual Velocity and Area
        """  
        Area_Actual = (Number * 0.25 * Hole_Diameter**2 * np.pi).to(ureg.inch**2)
        Velocity_Actual = (self.mdot/(self.rho * Area_Actual)).to(ureg.feet / ureg.second)
        self.Area_Actual = Area_Actual
        self.Velocity_Actual = Velocity_Actual
        return Velocity_Actual, Area_Actual           


# -------------- Design Inputs and Constants -------------- #    
#Constants
CD_drill = 0.7 #Constant for Sharp Edged Orifices
g0 = Q_(32.174, ureg.foot / ureg.second**2)
Prescott_pressure = Q_(12.04, ureg.force_pound / ureg.inch**2)
#Design Input
mdots = np.array([5.29, 2.21])/4 #LOX_CORE, FUEL_CORE
Film_Cooling = np.array([0.08, 0.08]) #Outer Film Cooling Percentage, Inner Film Cooling Percentage
di = 6.5 #Internal Diameter of Chamber
ri = di / 2 #Internal Radius of Chamber
Spacing = 0.75  #Spacing between center of impingement Holes
Rgamma_lox = 1.75  #Radial distance between centerline and LOX hole
Pressure_Drop_Fuel = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Drop_Lox = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Chamber = Q_(300, ureg.force_pound / ureg.inch**2) #Chamber Pressure Pretty Obvious
Doublet_Diameter_LOX = Q_(0.125, ureg.inch)  #Design choise for DOublet Diameter size (Need to look more properly into it as 1/4 holes might make vaporization time too long)\
Lox_Dewar_Pressure = Q_(22, ureg.force_pound / ureg.inch**2) + Prescott_pressure

# -------------- Prop Initialization -------------- #
LOX_Sat_Temp = LOX.Ts(p=Lox_Dewar_Pressure.magnitude )
LOX_Sat_Dens = LOX.ds(T=LOX_Sat_Temp)
OX_CORE = PROP(gamma=20., mdot=mdots[0], rho=LOX_Sat_Dens[0][0])
FUEL_CORE = PROP(gamma = 0, mdot = mdots[1], rho=51.15666) #gamma zero for this one because it's the initialized guess just making the FUEL CORE class requires it ( should change when moving to data classes)
OUT_FILM_C = PROP(gamma = 10, mdot = Film_Cooling[0]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
IN_FILM_C = PROP(gamma = -10, mdot = Film_Cooling[1]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)

# -------------- Function in DOublet to make my version of the SHitty Spike Contour -------------- #
Points = 1000 #also used in other plots if this section ends up getting deleted
x_profile,y_profile = spike_contour(Points)

# -------------- Code to Find Peaks (Largest diameter)  for Spike COntour (Will work with davids 2d code) -------------- #
Peaky, Peak_Point = max(y_profile), np.argmax(y_profile)
Peakx = x_profile[Peak_Point]

OX_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Lox)
FUEL_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
# -------------- Code to numericlaly find Fuel gamma based on the physical constraints and Ox gamma choice -------------- #
def func(gamma_FUEL):
    """Numerically solving function that has Momentum Balance on one side, and Tan Resultant angle in terms of gammas and physical Dimensions

    Args:
        gamma_FUEL (float): Angle off axial for fuel to impinge at

    Returns:
        float: Returns what the numerical solver goes for A.K.A zero
    """    
    return \
    -((ri + Peaky)/2 - Rgamma_lox -Spacing * np.sin(np.pi/2 + gamma_FUEL) / np.sin(np.radians(OX_CORE.gamma.magnitude) - gamma_FUEL) * np.cos(np.radians(90 - (OX_CORE.gamma.magnitude))))/ \
    (Peakx - Spacing * np.sin(np.pi/2 + gamma_FUEL) / np.sin(np.radians(OX_CORE.gamma.magnitude) - gamma_FUEL) * np.sin(np.radians(90 - (OX_CORE.gamma.magnitude)))) +\
    (OX_CORE.mdot.magnitude * OX_CORE.Velocity_init.magnitude * np.sin(np.deg2rad(OX_CORE.gamma.magnitude)) \
    + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_init.magnitude * np.sin(gamma_FUEL))/ \
    (OX_CORE.mdot.magnitude * OX_CORE.Velocity_init.magnitude * np.cos(np.deg2rad(OX_CORE.gamma.magnitude)) \
    + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_init.magnitude * np.cos(gamma_FUEL))    


gamma_FUEL_original = fsolve(func, np.deg2rad(5)) 
FUEL_CORE.gamma = Q_(np.rad2deg(gamma_FUEL_original[0]), ureg.degrees) 
print(f"Originally Solved Fuel Angle is {FUEL_CORE.gamma:.3f}")
FUEL_CORE.gamma = Q_(round(FUEL_CORE.gamma.magnitude * 2)/2, ureg.degrees) 
print(f"Newly Machinable Adjusted Fuel Angle is {FUEL_CORE.gamma:.3f}")


# -------------- Hole Size SHit -------------- #
OX_CORE.Number(Doublet_Diameter_LOX,CD_drill, Pressure_Chamber * (Pressure_Drop_Lox))
print(f'Number of Oxygen Doublet holes needed based on a {Doublet_Diameter_LOX} diameter is {OX_CORE.Hole_Number.magnitude} holes')
FUEL_CORE_Diameter = np.sqrt((FUEL_CORE.Area(CD_drill, Pressure_Chamber * (Pressure_Drop_Fuel)) / OX_CORE.Hole_Number.magnitude) * 4 / np.pi)
while True:
    Doublet_Diameter_Fuel , drill_size ,closest_index = drill_approximation(FUEL_CORE_Diameter.magnitude)
    Doublet_Diameter_Fuel = Q_(Doublet_Diameter_Fuel, ureg.inch)
    FUEL_CORE.Number((Doublet_Diameter_Fuel),CD_drill, Pressure_Chamber * (Pressure_Drop_Fuel))
    if FUEL_CORE.Hole_Number.magnitude == OX_CORE.Hole_Number.magnitude:
        break
    elif FUEL_CORE.Hole_Number.magnitude > OX_CORE.Hole_Number.magnitude:
        FUEL_CORE_Diameter += Q_(0.001,ureg.inch)
    else:
        FUEL_CORE_Diameter -= Q_(0.001,ureg.inch)
print(f"Closest drill size to {FUEL_CORE_Diameter :.5f~} is a diameter of {Doublet_Diameter_Fuel} with a drill size of {drill_size} .")
print(f'Number of Fuel Doublet holes needed based on a {Doublet_Diameter_Fuel} diameter is {FUEL_CORE.Hole_Number.magnitude} holes')
OX_CORE.Actual(Doublet_Diameter_LOX,OX_CORE.Hole_Number) #Initializing the actual function to use actual velocities
FUEL_CORE.Actual(Doublet_Diameter_Fuel,FUEL_CORE.Hole_Number) #Initializing the actual function to use actual velocities



# PLOTTING SHIT BELOW
# -------------- Constants and parameters -------------- #
gamma_lox = OX_CORE.gamma.magnitude  # degrees and making the variables work below since i made the matlab version of this first and converted to python with ChatGPT
gamma_fuel = FUEL_CORE.gamma.magnitude  # degrees
Chamber_Cowl_r = 0.5  # in
Past_Peak = 1.15 #some terrible constant for shitty chamber I made drawn


# -------------- Plotting the aerospike nozzle contour -------------- #
plt.figure()
plt.plot(x_profile, y_profile, linewidth=2) 
plt.axhline(0, linewidth=2, color="b")


# -------------- Calculations and Plotting for impingement point and Resultant Point -------------- #
x = Spacing * np.sin(np.radians(90 + gamma_fuel)) / np.sin(np.radians(gamma_lox - gamma_fuel)) * np.sin(np.radians(90 - gamma_lox)) 
y = Spacing * np.sin(np.radians(90 + gamma_fuel)) / np.sin(np.radians(gamma_lox - gamma_fuel)) * np.cos(np.radians(90 - gamma_lox)) + Rgamma_lox
plt.plot(x, y, "o")
yprime = (ri + Peaky) / 2
xprime = Peakx
plt.plot(xprime, yprime, "o")


# -------------- Plotting Horizontal Lines axial lines (References for angles) -------------- #
plt.axhline(Rgamma_lox, linestyle="--")
plt.axhline(Rgamma_lox + Spacing, linestyle="--")
plt.axhline(y, linestyle="--") #need y before this plot. THis is the horizontal line to represent the axial POV from the resultant angle


# -------------- Creating all linspaces needed for the following plots -------------- #
x_graph = np.linspace(0, max(x_profile), Points)
x_angled_lines = np.linspace(0, x, Points)  # Up to the impingement point for FUEL line
thetaRange4 = np.linspace(np.radians(90), 0, Points)


# -------------- Making and Plotting a Shitty Chamber Drawing  -------------- 
ChamberX = x_graph / x_graph[-1] * Peakx * Past_Peak
ChamberY = np.ones(len(ChamberX)) * ri
ChamberArcX = ChamberX[-1] + Chamber_Cowl_r * np.cos(thetaRange4)
ChamberArcY = ChamberY[-1] + Chamber_Cowl_r * (-1 + np.sin(thetaRange4))
Chamber_ContourX = np.concatenate([ChamberX, ChamberArcX])
Chamber_ContourY = np.concatenate([ChamberY, ChamberArcY])
plt.plot(Chamber_ContourX, Chamber_ContourY, "k", linewidth=2)


# -------------- Plotting the angled Prop Lines -------------- #
gamma_lox_line = np.tan(np.radians(gamma_lox)) * x_angled_lines + Rgamma_lox
plt.plot(x_angled_lines, gamma_lox_line,"g")
gamma_fuel_line = np.tan(np.radians(gamma_fuel)) * x_angled_lines + Rgamma_lox + Spacing
plt.plot(x_angled_lines, gamma_fuel_line,"r")

# -------------- Plotting the film cooling lines -------------- #
innercooling_line 
outercooling_line

gamma_lox_line = np.tan(np.radians(gamma_lox)) * x_angled_lines + Rgamma_lox
plt.plot(x_angled_lines, gamma_lox_line,"g")
gamma_fuel_line = np.tan(np.radians(gamma_fuel)) * x_angled_lines + Rgamma_lox + Spacing
plt.plot(x_angled_lines, gamma_fuel_line,"r")


# -------------- Plotting the resultant line -------------- #
tan_resultant =     (OX_CORE.mdot.magnitude * OX_CORE.Velocity_Actual.magnitude * np.sin(np.deg2rad(OX_CORE.gamma.magnitude)) \
    + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_Actual.magnitude * np.sin(np.deg2rad(FUEL_CORE.gamma.magnitude)))/ \
    (OX_CORE.mdot.magnitude * OX_CORE.Velocity_Actual.magnitude * np.cos(np.deg2rad(OX_CORE.gamma.magnitude)) \
    + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_Actual.magnitude * np.cos(np.deg2rad(FUEL_CORE.gamma.magnitude))) 
resultant_y_intercept = y - tan_resultant * x
ResultantX = x_graph / x_graph[-1] * (xprime - x) * Past_Peak + x
ResultantY = tan_resultant * ResultantX + resultant_y_intercept
plt.plot(ResultantX, ResultantY, "y", linewidth=2)


# -------------- Plotting the resultant error line And Error -------------- #
a,b,c = -tan_resultant,1,-resultant_y_intercept
Error = Q_(np.abs(a*xprime + b*yprime +c)/np.sqrt(a**2 + b**2),ureg.inch)
print(f'The following errors are gathered from rounding to nearest drill size, and half degree for angles:\
    \n\tFuel gamma was rounded from {np.rad2deg(gamma_FUEL_original[0]) :.4f} to {gamma_fuel}\
    \n\tLox Velocity was supposed to be {OX_CORE.Velocity_init :.4f~} but ended up at {OX_CORE.Velocity_Actual :.4f~} with the {Doublet_Diameter_LOX :.4f} holes\
    \n\tFuel Velocity was supposed to be {FUEL_CORE.Velocity_init :.4f~} but ended up at {FUEL_CORE.Velocity_Actual :.4f~} with the {Doublet_Diameter_Fuel :.4f} holes\
    \n\tWhich all resulted in missing the target by {Error :.4f~}')


# -------------- Plotting the arcs that show the angles -------------- #
arc_lox = patches.Arc((0, Rgamma_lox), x, x, 
                      angle=0, theta1=0, theta2=gamma_lox, color='green', label=f'Gamma_LOX: {gamma_lox} deg')
plt.gca().add_patch(arc_lox)
arc_fuel = patches.Arc((0, Rgamma_lox + Spacing), x, x, 
                       angle=0, theta1=gamma_fuel, theta2=0, color='red', label=f'Gamma_FUEL: {gamma_fuel} deg')
plt.gca().add_patch(arc_fuel)
delta = Q_((np.arctan(tan_resultant)) *180 / np.pi, ureg.degrees)
arc_delta = patches.Arc((x, y), 3*x, 3*x, 
                       angle=0, theta1=0, theta2=delta.magnitude, color='y', label=f'Resultant_Angle: {delta} deg')
plt.gca().add_patch(arc_delta)


# -------------- Plotting the angles onto the graph for easier readibility -------------- #
x_offset = -x/2  # Adjust as needed
y_offset_lox = -3*(y - Rgamma_lox)/4  
y_offset_fuel = 5*(Spacing + Rgamma_lox - y)/8  
x_offset_delta = 3*(xprime-x)/4  
y_offset_delta = 5*(yprime-y)/8
FUEL_CORE_gamma_writing = abs(FUEL_CORE.gamma.magnitude)
plt.text(x + x_offset,y + y_offset_lox, rf'$\gamma_{{\mathrm{{LOX}}}}: {OX_CORE.gamma.magnitude:.2f}^\circ$', color='green')
plt.text(x + x_offset,y + y_offset_fuel, rf'$\gamma_{{\mathrm{{FUEL}}}}: {FUEL_CORE_gamma_writing:.2f}^\circ$', color='red')
plt.text(x + x_offset_delta, y + y_offset_delta, rf'$\delta: {delta.magnitude:.2f}^\circ$', color='y')


# -------------- Extra Plotting Shit -------------- #
plt.legend(['Spike Contour', 'Centerline', 'Impingement Point', 'Aim Point', 'Gamma_(OX) Straight Line', 'Gamma_(FUEL) Straight Line','Resultant Straight Line',
             'Chamber Contour', f'Gamma_(OX) Angled Line {OX_CORE.gamma :.3f~}', 
             f'Gamma_(FUEL) Angled Line {FUEL_CORE.gamma :.3f~}', f'Resultant Line {delta :.3f~}'], loc="best")
plt.xlabel('Distance Along Engine Axis (inches)')
plt.ylabel('Radius (inches)')
plt.axis('equal')
plt.title('Side View Contour of an Aerospike Nozzle')
plt.grid(True)
plt.show()


# injector_cad_write()