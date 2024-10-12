from matplotlib import pyplot as plt, patches
from Doublet_Functions import spike_contour
from InjectorCad import injector_cad_write
from Drill import drill_approximation
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from fluid import Q_, ureg, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox, PROPFLOWS
import fluid
import numpy as np
import os
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
# -------------- Design Inputs -------------- #    
di = 6.5 #Internal Diameter of Chamber
ri = di / 2 #Internal Radius of Chamber
Spacing = 0.55  #Spacing between center of impingement Holes
Rgamma_lox = 1.65  #Radial distance between centerline and LOX hole
mdots = np.array([5.29, 2.21])/4 #LOX_CORE, FUEL_CORE lbm/s
Film_Cooling = np.array([0.05, 0.05]) #Outer Film Cooling Percentage, Inner Film Cooling Percentage
FilmCoolingSpacing = np.array([.25, .25]) #inches Inner, Outer
Pressure_Chamber = Q_(300, ureg.force_pound / ureg.inch**2) #Chamber Pressure Pretty Obvious
Doublet_Diameter_LOX = Q_(0.0625, ureg.inch)  #Design choise for DOublet Diameter size (Need to look more properly into it as 1/4 holes might make vaporization time too long)\
gammas = np.array([25.,0,30,-30]) #Lox, fuel, outer FC, inner FC  #All Angles from axial
Lox_Dewar_Pressure = 22

# -------------- Prop Initialization -------------- #
OX_CORE,FUEL_CORE,OUT_FILM_C,IN_FILM_C = PROPFLOWS(mdots,Film_Cooling,gammas,Lox_Dewar_Pressure)

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
        FUEL_CORE_Diameter += Q_(0.0005,ureg.inch)
    else:
        FUEL_CORE_Diameter -= Q_(0.0005,ureg.inch)
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


# -------------- Creating all linspaces needed for the following plots -------------- #
x_graph = np.linspace(0, max(x_profile), Points)
x_angled_lines = np.linspace(0, x, Points)  # Up to the impingement point for FUEL line
print(FilmCoolingSpacing[0]* np.tan( np.pi / 2  - np.radians(IN_FILM_C.gamma.magnitude)))
x_filmcooling_inner = np.linspace(0, FilmCoolingSpacing[0]* np.tan( np.pi / 2  - np.abs(np.radians(IN_FILM_C.gamma.magnitude))), Points)
x_filmcooling_outer = np.linspace(0, FilmCoolingSpacing[1]* np.tan( np.pi / 2  - np.abs(np.radians(OUT_FILM_C.gamma.magnitude))), Points)
thetaRange4 = np.linspace(np.radians(90), 0, Points)


# -------------- Making and Plotting a Shitty Chamber Drawing  -------------- 
ChamberX = x_graph / x_graph[-1] * Peakx * Past_Peak
ChamberY = np.ones(len(ChamberX)) * ri
ChamberArcX = ChamberX[-1] + Chamber_Cowl_r * np.cos(thetaRange4)
ChamberArcY = ChamberY[-1] + Chamber_Cowl_r * (-1 + np.sin(thetaRange4))
Chamber_ContourX = np.concatenate([ChamberX, ChamberArcX])
Chamber_ContourY = np.concatenate([ChamberY, ChamberArcY])
plt.plot(Chamber_ContourX, Chamber_ContourY, "k", linewidth=2)

# -------------- Plotting Horizontal Lines axial lines (References for angles) -------------- #
plt.axhline(Rgamma_lox, linestyle="--")
plt.axhline(Rgamma_lox + Spacing, linestyle="--")
plt.axhline(y, linestyle="--") #need y before this plot. THis is the horizontal line to represent the axial POV from the resultant angle
plt.axhline((y_profile[0] + FilmCoolingSpacing[0]) , linestyle="dotted")
plt.axhline((ChamberY[0] - FilmCoolingSpacing[1]) , linestyle="dotted")

# -------------- Plotting the angled Prop Lines -------------- #
gamma_lox_line = np.tan(np.radians(gamma_lox)) * x_angled_lines + Rgamma_lox
plt.plot(x_angled_lines, gamma_lox_line,"g")
gamma_fuel_line = np.tan(np.radians(gamma_fuel)) * x_angled_lines + Rgamma_lox + Spacing
plt.plot(x_angled_lines, gamma_fuel_line,"r")

# -------------- Plotting the film cooling lines -------------- #
innercooling_line = np.tan(np.radians(IN_FILM_C.gamma.magnitude)) * x_filmcooling_inner + y_profile[0] + FilmCoolingSpacing[0] 
plt.plot(x_filmcooling_inner, innercooling_line,"r")
outercooling_line = np.tan(np.radians(OUT_FILM_C.gamma.magnitude)) * x_filmcooling_outer + ChamberY[0] - FilmCoolingSpacing[1] 
plt.plot(x_filmcooling_outer, outercooling_line,"r")


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
arc_filmcooling_inner = patches.Arc((0, (y_profile[0] + FilmCoolingSpacing[0])), 1.5*x, 1.5*x, 
                      angle=0, theta1=IN_FILM_C.gamma.magnitude, theta2=0, color='red', label=f'Gamma_Filmcooling_Inner: {IN_FILM_C.gamma.magnitude} deg')
plt.gca().add_patch(arc_filmcooling_inner)
arc_filmcooling_outer = patches.Arc((0, (ChamberY[0] - FilmCoolingSpacing[1])), 1.5*x, 1.5*x, 
                      angle=0, theta1=0, theta2=OUT_FILM_C.gamma.magnitude, color='red', label=f'Gamma_FilmCooling_Outer: {OUT_FILM_C.gamma.magnitude} deg')
plt.gca().add_patch(arc_filmcooling_outer)


# -------------- Plotting the angles onto the graph for easier readibility -------------- #
x_offset = -x/2  # Adjust as needed
y_offset_lox = -3*(y - Rgamma_lox)/4  
y_offset_fuel = 5*(Spacing + Rgamma_lox - y)/8  
x_offset_delta = 3*(xprime-x)/4  
y_offset_delta = 5*(yprime-y)/8
x_offset_Filmcooling = x
y_offset_Filmcooling_Inner = -2*FilmCoolingSpacing[0]/4
y_offset_Filmcooling_Outer =  FilmCoolingSpacing[1]/4
FUEL_CORE_gamma_writing = abs(FUEL_CORE.gamma.magnitude)
plt.text(x + x_offset,y + y_offset_lox, rf'$\gamma_{{\mathrm{{OX}}}}: {OX_CORE.gamma.magnitude:.2f}^\circ$', color='green')
plt.text(x + x_offset,y + y_offset_fuel, rf'$\gamma_{{\mathrm{{F}}}}: {FUEL_CORE_gamma_writing:.2f}^\circ$', color='red')
plt.text(x + x_offset_delta, y + y_offset_delta, rf'$\delta: {delta.magnitude:.2f}^\circ$', color='y')
plt.text(x_offset_Filmcooling, (y_profile[0] + FilmCoolingSpacing[0]) + y_offset_Filmcooling_Inner, rf'$\gamma_{{\mathrm{{FCI}}}}: {IN_FILM_C.gamma.magnitude:.2f}^\circ$', color='red')
plt.text(x_offset_Filmcooling, (ChamberY[0] - FilmCoolingSpacing[1]) + y_offset_Filmcooling_Outer, rf'$\gamma_{{\mathrm{{FCO}}}}: {OUT_FILM_C.gamma.magnitude:.2f}^\circ$', color='red')


# -------------- Extra Plotting Shit -------------- #
plt.legend(['Spike Contour', 'Centerline', 'Impingement Point', 'Aim Point','Chamber Contour', 'Gamma_(OX) Straight Line', 'Gamma_(FUEL) Straight Line','Resultant Straight Line', 
            '(Film Cooling Inner) Straight Line', 'Film Cooling Outer Straight Line',
              f'Gamma_(OX) Angled Line {OX_CORE.gamma :.3f~}', 
             f'Gamma_(FUEL) Angled Line {FUEL_CORE.gamma :.3f~}', 
             f'Gamma_FilmCooling_Inner {IN_FILM_C.gamma :.3f~}',
             f'Gamma_FilmCooling_Outer {OUT_FILM_C.gamma :.3f~}',
             f'Resultant Line {delta :.3f~}'], loc="upper right", bbox_to_anchor=(1.05,.85), mode="expand")
plt.subplots_adjust(left=0.06)
plt.subplots_adjust(bottom=0.069)
plt.subplots_adjust(right=0.7420)
plt.subplots_adjust(top=0.9)
plt.xlabel('Distance Along Engine Axis (inches)')
plt.ylabel('Radius (inches)')
plt.axis('equal')
plt.title('Side View Contour of an Aerospike Nozzle')
plt.grid(True)
plt.show()


# injector_cad_write()