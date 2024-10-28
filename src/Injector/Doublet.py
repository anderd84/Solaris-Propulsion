from matplotlib import pyplot as plt, patches
from Injector.Doublet_Functions import spike_contour
from Injector.InjectorCad import injector_cad_write
from Injector.Drill import drill_approximation
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from fluids.fluid import Q_, ureg, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox, PROPFLOWS, Chamber_Press
import numpy as np
from icecream import ic
from texttable import Texttable

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
di = 5.5 #Internal Diameter of Chamber
ri = di / 2 #Internal Radius of Chamber
Spacing = 0.55  #Spacing between centear of impingement Holes
Rgamma_lox = 1.25  #Radial distance between centerline and LOX hole
Film_Cooling = np.array([0.15, 0.00]) #Outer Film Cooling Percentage, Inner Film Cooling Percentage
FilmCoolingSpacing = np.array([.60, .25]) #inches Inner, Outer
Pressure_Chamber = Q_(Chamber_Press, ureg.force_pound / ureg.inch**2) #Chamber Pressure Pretty Obvious
Doublet_Diameter_LOX = Q_(0.0625, ureg.inch)  #Design choise for DOublet Diameter size (Need to look more properly into it as 1/4 holes might make vaporization time too long)\
gammas = np.array([25.,0,25,-25]) #Lox, fuel, outer FC, inner FC  #All Angles from axial
Lox_Dewar_Pressure = 22
AirTemperature = Q_(70, ureg.degF)
AirPressure = Q_(12.2, ureg.force_pound / ureg.inch**2)
FuelName = "RP-1" 
Points = 1000 #also used in other plots if this section ends up getting deleted

def InjectorParameters(ri: float ,Spacing: float , Rgamma_lox: float, Film_Cooling: float, FilmCoolingSpacing: float, Pressure_Chamber: float,Doublet_Diameter_LOX: float
                       , gammas: float,Lox_Dewar_Pressure: float, AirTemperature: float, AirPressure: float ,FuelName, spike_contour: float, Points: float):
    # -------------- Prop Initialization -------------- #
    OX_CORE,FUEL_CORE,OUT_FILM_C,IN_FILM_C,viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f = PROPFLOWS(Film_Cooling,gammas,Lox_Dewar_Pressure, AirTemperature, AirPressure ,  FuelName)

    # -------------- Function in DOublet to make my version of the SHitty Spike Contour -------------- #
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
    gamma_FUEL_original = Q_(np.rad2deg(gamma_FUEL_original[0]), ureg.degrees) 
#    FUEL_CORE.gamma = Q_(np.rad2deg(gamma_FUEL_original.magnitude), ureg.degrees) 
    FUEL_CORE.gamma = Q_(round(gamma_FUEL_original.magnitude * 2)/2, ureg.degrees) 


    # -------------- Main Doublets Holes -------------- #
    # Initialize variables to track the closest match
    OX_CORE.Number(Doublet_Diameter_LOX, CD_drill, Pressure_Chamber * Pressure_Drop_Lox)
    FUEL_CORE_Diameter = np.sqrt((FUEL_CORE.Area(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel) / OX_CORE.Hole_Number.magnitude) * 4 / np.pi)
    target_hole_number = OX_CORE.Hole_Number.magnitude  # Target number of holes
    # Set initial values for the closest diameter match
    closest_diameter = FUEL_CORE_Diameter
    FUEL_CORE.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    closest_hole_diff = 10
    iterated = 100
    while iterated > 0:
        # Calculate new diameter and corresponding hole count
        Doublet_Diameter_Fuel, drill_size, closest_index = drill_approximation(FUEL_CORE_Diameter.magnitude)
        Doublet_Diameter_Fuel = Q_(Doublet_Diameter_Fuel, ureg.inch)
        FUEL_CORE.Number(Doublet_Diameter_Fuel, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
        # Calculate current difference in hole numbers
        hole_diff = np.abs(FUEL_CORE.Hole_Number.magnitude - target_hole_number)
        # Update the closest match if this iteration is better
        if hole_diff < closest_hole_diff:
            closest_diameter = Doublet_Diameter_Fuel
            closest_hole_diff = hole_diff
        # Break if exact match is found
        if hole_diff == 0:
            break
        # Adjust diameter based on hole count difference
        if FUEL_CORE.Hole_Number.magnitude > target_hole_number:
            FUEL_CORE_Diameter += Q_(0.0005, ureg.inch)
        else:
            FUEL_CORE_Diameter -= Q_(0.0005, ureg.inch)
        iterated -= 1
    # After exiting loop, set FUEL_CORE to the closest diameter found
    # Initialize the actual values for the selected diameter
    FUEL_CORE.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    OX_CORE.Actual(Doublet_Diameter_LOX, OX_CORE.Hole_Number)
    FUEL_CORE.Actual(closest_diameter, OX_CORE.Hole_Number)
    fuel_Doublet = [FUEL_CORE_Diameter, closest_diameter, drill_size]



    # -------------- Outer Film Cooling Holes -------------- #
    RFilmCooling = ri - FilmCoolingSpacing[0]
    OuterFC_Number_Guess = np.ceil(2 * np.pi * RFilmCooling / 0.3)  # NASA spec: 0.3 inches between holes
   # Initial estimate for Outer Film Cooling hole diameter
    OuterFC_Diameter = np.sqrt((OUT_FILM_C.Area(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel) / OuterFC_Number_Guess) * 4 / np.pi)
    target_hole_number = OuterFC_Number_Guess  # Target number of holes for the outer film cooling
    # Set initial closest diameter and hole difference
    closest_diameter = OuterFC_Diameter
    OUT_FILM_C.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    closest_hole_diff = 10
    iterated = 100
    while iterated > 0:
        # Determine the diameter and corresponding hole count
        Doublet_Diameter_OuterFC, drill_size, closest_index = drill_approximation(OuterFC_Diameter.magnitude)
        Doublet_Diameter_OuterFC = Q_(Doublet_Diameter_OuterFC, ureg.inch)
        OUT_FILM_C.Number(Doublet_Diameter_OuterFC, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)

        # Calculate the difference in hole numbers
        hole_diff = np.abs(OUT_FILM_C.Hole_Number.magnitude - target_hole_number)
        # Update closest match if current iteration is better
        if hole_diff < closest_hole_diff:
            closest_diameter = Doublet_Diameter_OuterFC
            closest_hole_diff = hole_diff
        # Break if exact match is found
        if hole_diff == 0:
            break
        # Adjust diameter based on hole count difference
        if OUT_FILM_C.Hole_Number.magnitude > target_hole_number:
            OuterFC_Diameter += Q_(0.0001, ureg.inch)
        else:
            OuterFC_Diameter -= Q_(0.0001, ureg.inch)
        iterated -= 1

    # After loop, set the closest diameter found
    OUT_FILM_C.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    # Initialize the actual values for the selected diameter
    OUT_FILM_C.Actual(closest_diameter, OUT_FILM_C.Hole_Number)
    film_Cool_Doublet = [OuterFC_Diameter, closest_diameter, drill_size]

    # DROPLET SIZING
    # R. A. DlCKERSON method
    D_f_Dickerson = Q_( 1e5 * np.power(Doublet_Diameter_Fuel.magnitude, 0.27) *  np.power(Doublet_Diameter_LOX.magnitude, 0.023) / \
        (  np.power(FUEL_CORE.Velocity_Actual.magnitude, 0.74) * np.power(OX_CORE.Velocity_Actual.magnitude, 0.33)), ureg.micron   )
    P_dynamic_OX = (0.5 * OX_CORE.rho * OX_CORE.Velocity_Actual**2).to(ureg.psi)
    P_dynamic_FUEL = (0.5 * FUEL_CORE.rho * FUEL_CORE.Velocity_Actual**2).to(ureg.psi)
    # Now use these in the NASA SP-8089 droplet size calculation
    P_D = (P_dynamic_OX / P_dynamic_FUEL).magnitude
    viscosity_shellwax = Q_(2.69e-3, ureg.pound / (ureg.foot * ureg.second))  # lbm/(ft-sec)
    rho_shellwax = Q_(47.7, ureg.pound / ureg.foot**3)  # lbm/ft^3
    SurfaceTens_shellwax = Q_(17, ureg.dyne / ureg.centimeter)  # dynes/cm
    viscosity_f = viscosity_f.to(ureg.pound / (ureg.foot * ureg.second))
    SurfaceTens_f = Q_(SurfaceTens_f.to(ureg.dyne / ureg.centimeter))
    # Assuming `FUEL_CORE.Viscosity`, `FUEL_CORE.Density`, and `FUEL_CORE.Surface_Tension` are defined
    # and contain the viscosity, density, and surface tension of the fuel, respectively.
    Pc_Pj = 2 #Ratio of the Core Max Velocity / Average Velocity Roughly 2
    K_prop = np.power(
        ((viscosity_f * SurfaceTens_f /FUEL_CORE.rho) / (viscosity_shellwax * SurfaceTens_shellwax /rho_shellwax)).magnitude, 1 / 4 )
    # Final NASA SP-8089 droplet size calculation
    D_f_NASA = Q_(
        2.9e4 * np.power(FUEL_CORE.Velocity_Actual.magnitude, -0.766) * np.power(Pc_Pj, -0.65) * 
        np.power(P_D, 0.165) * np.power(Doublet_Diameter_Fuel.magnitude, 0.293) *
        np.power(Doublet_Diameter_LOX.magnitude / Doublet_Diameter_Fuel.magnitude, 0.023) * K_prop,
        ureg.micron)
    # VAPORIZATION TIME & CHAMBER LENGTH
    B = 10 # B for most hydrocarbon fuels are between 5 & 20.
    Vaporize_time_Dickerson = (FUEL_CORE.rho * specific_heat_p_f * (D_f_Dickerson/2)**2 / (2* thermal_conductivity_f * np.log(1 + B))).to(ureg.second)
    Vaporize_time_NASA = (FUEL_CORE.rho * specific_heat_p_f * (D_f_NASA/2)**2 / (2* thermal_conductivity_f * np.log(1 + B))).to(ureg.second)
    gamma_lox_rad = np.deg2rad(OX_CORE.gamma.magnitude)
    gamma_fuel_rad = np.deg2rad(FUEL_CORE.gamma.magnitude)
    # Calculate momentum components
    momentum_x = (OX_CORE.mdot * OX_CORE.Velocity_Actual * np.cos(gamma_lox_rad) +
                  FUEL_CORE.mdot * FUEL_CORE.Velocity_Actual * np.cos(gamma_fuel_rad))
    momentum_y = (OX_CORE.mdot * OX_CORE.Velocity_Actual * np.sin(gamma_lox_rad) +
                  FUEL_CORE.mdot * FUEL_CORE.Velocity_Actual * np.sin(gamma_fuel_rad))
    Resultant_angle = np.arctan2(momentum_y,momentum_x)
    # Calculate total resultant impingement velocity from x and y components
    V_impingement = (np.sqrt(momentum_x**2 + momentum_y**2) / (OX_CORE.mdot + FUEL_CORE.mdot)).to(ureg.feet / ureg.second)
    Travel_Length_Dickerson = (V_impingement *Vaporize_time_Dickerson).to(ureg.inch)
    Travel_Length_NASA = (V_impingement *Vaporize_time_NASA).to(ureg.inch)
    Chamber_Length_Dickerson = Travel_Length_Dickerson * np.cos(Resultant_angle)
    Chamber_Length_NASA = Travel_Length_NASA * np.cos(Resultant_angle)
    Dickerson = [D_f_Dickerson, Vaporize_time_Dickerson, Travel_Length_Dickerson, Chamber_Length_Dickerson]
    NASA = [D_f_NASA, Vaporize_time_NASA, Travel_Length_NASA, Chamber_Length_NASA]
 
    return OX_CORE, FUEL_CORE, IN_FILM_C, OUT_FILM_C, Dickerson, NASA, fuel_Doublet, film_Cool_Doublet, Doublet_Diameter_Fuel, gamma_FUEL_original, x_profile, y_profile, Peakx, Peaky


OX_CORE, FUEL_CORE, IN_FILM_C, OUT_FILM_C, Dickerson, NASA, fuel_Doublet, film_Cool_Doublet, Doublet_Diameter_Fuel, gamma_FUEL_original, x_profile, y_profile, Peakx, Peaky = InjectorParameters(ri ,Spacing, Rgamma_lox, Film_Cooling, FilmCoolingSpacing, Pressure_Chamber,Doublet_Diameter_LOX, gammas,Lox_Dewar_Pressure, AirTemperature, AirPressure ,FuelName, spike_contour,Points)
# PLOTTING SHIT BELOW



fuel_Core_Diameter, closest_Fuel_Diameter, fuel_Doublet_Drill = fuel_Doublet
film_Cool_CORE_Diameter, closest_Film_Cool_Diameter, film_Cool_Doublet_Drill = film_Cool_Doublet
[D_f_Dickerson, Vaporize_time_Dickerson, Travel_Length_Dickerson, Chamber_Length_Dickerson] = Dickerson
[D_f_NASA, Vaporize_time_NASA, Travel_Length_NASA, Chamber_Length_NASA] = NASA
_, LOX_doublet_drill_size, _ = drill_approximation(Doublet_Diameter_LOX.magnitude)
angle_off_axial_fuel_original = gamma_FUEL_original.magnitude  # Strip units for display
angle_off_axial_ox = OX_CORE.gamma.magnitude  # Strip units for display
angle_off_axial_fuel = FUEL_CORE.gamma.magnitude  # Strip units for display
angle_off_axial_film = OUT_FILM_C.gamma.magnitude  # Strip units for display


Injector_Parameters = Texttable()
Injector_Parameters.add_rows([["",'Oxidizer', 'Fuel', 'Film Cooling'], ["Prop Type","LOX", FuelName, FuelName]])
Injector_Parameters.add_row(["Angle off Axial (Deg)", angle_off_axial_ox, angle_off_axial_fuel, angle_off_axial_film])
Injector_Parameters.add_row(["Doublet Size", Doublet_Diameter_LOX,closest_Fuel_Diameter, closest_Film_Cool_Diameter ])
Injector_Parameters.add_row(["Drill Size", LOX_doublet_drill_size,fuel_Doublet_Drill, film_Cool_Doublet_Drill ])
Injector_Parameters.add_row(["Hole Number",OX_CORE.Number_Actual , FUEL_CORE.Number_Actual, OUT_FILM_C.Number_Actual ])
rounded_Table = Texttable()
rounded_Table.add_rows([["",'Unrounded Values', 'Rounded Values'], ["Fuel Angle off Axial (Deg)", angle_off_axial_fuel_original, angle_off_axial_fuel]])
rounded_Table.add_row(["Oxidizer Velocity", f"{OX_CORE.Velocity_init:.3f~}", f"{OX_CORE.Velocity_Actual:.3f~} "],)
rounded_Table.add_row(["Fuel Velocity", f"{FUEL_CORE.Velocity_init:.3f~}", f"{FUEL_CORE.Velocity_Actual:.3f~} "])


Vaporization = Texttable()
Vaporization.add_rows([["",'Dickerson Method', 'NASA Method'], ["Droplet Size (microns)", f"{D_f_Dickerson:.3f~} ", f"{D_f_NASA:.3f~} "]])
Vaporization.add_row(["Required Vaporization Time", f"{Vaporize_time_Dickerson:.3f~} ", f"{Vaporize_time_NASA:.3f~} "],)
Vaporization.add_row(["Required Travel Length", f"{Travel_Length_Dickerson:.3f~} ", f"{Travel_Length_NASA:.3f~} "])
Vaporization.add_row(["Required Chamber Length", f"{Chamber_Length_Dickerson:.3f~} ", f"{Chamber_Length_NASA:.3f~} "])

print(Injector_Parameters.draw())
print(rounded_Table.draw())
print(Vaporization.draw())

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
#print(FilmCoolingSpacing[0]* np.tan( np.pi / 2  - np.radians(IN_FILM_C.gamma.magnitude)))
x_filmcooling_inner = np.linspace(0, FilmCoolingSpacing[1]* np.tan( np.pi / 2  - np.abs(np.radians(IN_FILM_C.gamma.magnitude))), Points)
x_filmcooling_outer = np.linspace(0, FilmCoolingSpacing[0]* np.tan( np.pi / 2  - np.abs(np.radians(OUT_FILM_C.gamma.magnitude))), Points)
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
plt.axhline((y_profile[0] + FilmCoolingSpacing[1]) , linestyle="dotted")
plt.axhline((ri - FilmCoolingSpacing[0]) , linestyle="dotted")

# -------------- Plotting the angled Prop Lines -------------- #
gamma_lox_line = np.tan(np.radians(gamma_lox)) * x_angled_lines + Rgamma_lox
plt.plot(x_angled_lines, gamma_lox_line,"g")
gamma_fuel_line = np.tan(np.radians(gamma_fuel)) * x_angled_lines + Rgamma_lox + Spacing
plt.plot(x_angled_lines, gamma_fuel_line,"r")

# -------------- Plotting the film cooling lines -------------- #
innercooling_line = np.tan(np.radians(IN_FILM_C.gamma.magnitude)) * x_filmcooling_inner + y_profile[0] + FilmCoolingSpacing[1] 
plt.plot(x_filmcooling_inner, innercooling_line,"r")
outercooling_line = np.tan(np.radians(OUT_FILM_C.gamma.magnitude)) * x_filmcooling_outer + ri - FilmCoolingSpacing[0] 
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
print(f'The errors from rounding to nearest drill size, and half degree for angles: resulted in missing the target by {Error :.4f~}')


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
arc_filmcooling_inner = patches.Arc((0, (y_profile[0] + FilmCoolingSpacing[1])), 1.5*x, 1.5*x, 
                      angle=0, theta1=IN_FILM_C.gamma.magnitude, theta2=0, color='red', label=f'Gamma_Filmcooling_Inner: {IN_FILM_C.gamma.magnitude} deg')
plt.gca().add_patch(arc_filmcooling_inner)
arc_filmcooling_outer = patches.Arc((0, (ri - FilmCoolingSpacing[0])), 1.5*x, 1.5*x, 
                      angle=0, theta1=0, theta2=OUT_FILM_C.gamma.magnitude, color='red', label=f'Gamma_FilmCooling_Outer: {OUT_FILM_C.gamma.magnitude} deg')
plt.gca().add_patch(arc_filmcooling_outer)


# -------------- Plotting the angles onto the graph for easier readibility -------------- #
x_offset = -x/2  # Adjust as needed
y_offset_lox = -3*(y - Rgamma_lox)/4  
y_offset_fuel = 5*(Spacing + Rgamma_lox - y)/8  
x_offset_delta = 3*(xprime-x)/4  
y_offset_delta = 5*(yprime-y)/8
x_offset_Filmcooling = x
y_offset_Filmcooling_Inner = -2*FilmCoolingSpacing[1]/4
y_offset_Filmcooling_Outer =  FilmCoolingSpacing[0]/4
FUEL_CORE_gamma_writing = abs(FUEL_CORE.gamma.magnitude)
plt.text(x + x_offset,y + y_offset_lox, rf'$\gamma_{{\mathrm{{OX}}}}: {OX_CORE.gamma.magnitude:.2f}^\circ$', color='green')
plt.text(x + x_offset,y + y_offset_fuel, rf'$\gamma_{{\mathrm{{F}}}}: {FUEL_CORE_gamma_writing:.2f}^\circ$', color='red')
plt.text(x + x_offset_delta, y + y_offset_delta, rf'$\delta: {delta.magnitude:.2f}^\circ$', color='y')
plt.text(x_offset_Filmcooling, (y_profile[0] + FilmCoolingSpacing[1]) + y_offset_Filmcooling_Inner, rf'$\gamma_{{\mathrm{{FCI}}}}: {IN_FILM_C.gamma.magnitude:.2f}^\circ$', color='red')
plt.text(x_offset_Filmcooling, (ri - FilmCoolingSpacing[0]) + y_offset_Filmcooling_Outer, rf'$\gamma_{{\mathrm{{FCO}}}}: {OUT_FILM_C.gamma.magnitude:.2f}^\circ$', color='red')


# -------------- Extra Plotting Shit -------------- #
plt.legend(['Spike Contour', 'Centerline', 'Impingement Point', 'Aim Point','Chamber Contour', 'Gamma_(OX) Straight Line', 'Gamma_(FUEL) Straight Line',
            'Resultant Straight Line', '(Film Cooling Inner) Straight Line', 'Film Cooling Outer Straight Line',
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