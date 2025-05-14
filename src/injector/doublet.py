from matplotlib import pyplot as plt, patches
from scipy.optimize import fsolve
from fluids.fluid import  PROPFLOWS
import numpy as np
from icecream import ic
from texttable import Texttable

from injector.drill import drill_approximation
from general.units import Q_, unitReg
import general.design as DESIGN

OF = DESIGN.OFratio
TotalMdot = DESIGN.totalmdot
Chamber_Press = DESIGN.chamberPressure

def initialize_prop_flows(Film_Cooling, oxImpingeAngle, fuelInitalImpingeAngle, filmImpingeAngle,
                          Lox_Dewar_Pressure, CD_drill,Pressure_Drop_Lox, Pressure_Drop_Fuel, AirTemperature, AirPressure, FuelName, Pressure_Chamber: float):
    """Initialize oxidizer and fuel core properties."""

    OX_CORE,FUEL_CORE,OUT_FILM_C,viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f = PROPFLOWS(
        Film_Cooling,oxImpingeAngle, fuelInitalImpingeAngle, filmImpingeAngle,Lox_Dewar_Pressure, AirTemperature, AirPressure ,  FuelName)
    OX_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Lox)
    FUEL_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    return OX_CORE,FUEL_CORE,OUT_FILM_C,viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f


# Calculation of fuel Impingement angle
def calculate_fuel_impinge_angle(OX_CORE, FUEL_CORE, AimY, Rgamma_lox, Spacing, AimX):
    """Calculate the fuel impinge angle based on physical constraints."""
    # -------------- Code to numericlaly find Fuel gamma based on the physical constraints and Ox gamma choice -------------- #
    def func(gamma_FUEL):
        """Numerically solving function that has Momentum Balance on one side, and Tan Resultant angle in terms of gammas and physical Dimensions

        Args:
            gamma_FUEL (float): Angle off axial for fuel to impinge at

        Returns:
            float: Returns what the numerical solver goes for A.K.A zero
        """    
        Num  = -(AimY - Rgamma_lox.magnitude -Spacing.magnitude * np.sin(np.pi/2 + gamma_FUEL) / np.sin(np.radians(OX_CORE.gamma.magnitude) - gamma_FUEL) * np.cos(np.radians(90 - (OX_CORE.gamma.magnitude))))
        Den1 = (AimX - Spacing.magnitude * np.sin(np.pi/2 + gamma_FUEL) / np.sin(np.radians(OX_CORE.gamma.magnitude) - gamma_FUEL) * np.sin(np.radians(90 - (OX_CORE.gamma.magnitude)))) 
        Num2 = (OX_CORE.mdot.magnitude * OX_CORE.Velocity_init.magnitude * np.sin(np.deg2rad(OX_CORE.gamma.magnitude))) 
        Num3 = (FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_init.magnitude * np.sin(gamma_FUEL))
        Den2 = (OX_CORE.mdot.magnitude * OX_CORE.Velocity_init.magnitude * np.cos(np.deg2rad(OX_CORE.gamma.magnitude)))
        Den3 = (FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_init.magnitude * np.cos(gamma_FUEL))
        
        return Num/(Den1)  + (Num2 + Num3)/(Den2 + Den3)
    
    gamma_FUEL_original = fsolve(func, np.deg2rad(5)) 
    gamma_FUEL_original = Q_(np.rad2deg(gamma_FUEL_original[0]), unitReg.degrees) 
    FUEL_CORE.gamma = Q_(round(gamma_FUEL_original.magnitude * 2)/2, unitReg.degrees) 
    return gamma_FUEL_original, FUEL_CORE


def calculate_fuel_diameters(OX_CORE, FUEL_CORE, Pressure_Chamber, Doublet_Diameter_LOX, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox):
    # -------------- Main Doublets Holes -------------- #
    # Initialize variables to track the closest match
    OX_CORE.Number(Doublet_Diameter_LOX, CD_drill, Pressure_Chamber * Pressure_Drop_Lox)
    FUEL_CORE_Diameter = np.sqrt((FUEL_CORE.Area(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel) / OX_CORE.Hole_Number.magnitude) * 4 / np.pi)
    Original_Fuel_Diameter = FUEL_CORE_Diameter
    target_hole_number = OX_CORE.Hole_Number.magnitude  # Target number of holes
    
    # Set initial values for the closest diameter match
    closest_diameter = FUEL_CORE_Diameter
    FUEL_CORE.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    closest_hole_diff = 10
    iterated = 100

    while iterated > 0:
        # Calculate new diameter and corresponding hole count
        Doublet_Diameter_Fuel, drill_size, closest_index = drill_approximation(FUEL_CORE_Diameter.magnitude)
        Doublet_Diameter_Fuel = Q_(Doublet_Diameter_Fuel, unitReg.inch)
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
            FUEL_CORE_Diameter += Q_(0.0005, unitReg.inch)
        else:
            FUEL_CORE_Diameter -= Q_(0.0005, unitReg.inch)
        iterated -= 1
    
    # After exiting loop, set FUEL_CORE to the closest diameter found
    # Initialize the actual values for the selected diameter
    FUEL_CORE.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    OX_CORE.Actual(Doublet_Diameter_LOX, OX_CORE.Hole_Number)
    FUEL_CORE.Actual(closest_diameter, OX_CORE.Hole_Number)
    fuel_Doublet = [ Original_Fuel_Diameter, closest_diameter, drill_size]
    return fuel_Doublet


def calculate_film_cooling_diameters(OUT_FILM_C, FilmCoolingSpacing, Pressure_Chamber, CD_drill, Pressure_Drop_Fuel, chamberatInjectorRadius, manual_override):



    # -------------- Outer Film Cooling Holes -------------- #
    RFilmCooling = chamberatInjectorRadius - FilmCoolingSpacing 
    ic(RFilmCooling)
    OuterFC_Number_Guess = np.ceil(2 * np.pi * RFilmCooling / Q_(0.35, unitReg.inch))  # NASA spec: 0.3 inches between holes
   # Initial estimate for Outer Film Cooling hole diameter
    OuterFC_Diameter = np.sqrt((OUT_FILM_C.Area(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel) / OuterFC_Number_Guess) * 4 / np.pi)
    Original_Film_Cooling_Diameter = OuterFC_Diameter
    target_hole_number = OuterFC_Number_Guess  # Target number of holes for the outer film cooling
    # Set initial closest diameter and hole difference
    closest_diameter = OuterFC_Diameter
    OUT_FILM_C.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    closest_hole_diff = 10
    iterated = 100
    while iterated > 0:
        # Determine the diameter and corresponding hole count
        Doublet_Diameter_OuterFC, drill_size, closest_index = drill_approximation(OuterFC_Diameter.magnitude)
        Doublet_Diameter_OuterFC = Q_(Doublet_Diameter_OuterFC, unitReg.inch)
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
            OuterFC_Diameter += Q_(0.0001, unitReg.inch)
        else:
            OuterFC_Diameter -= Q_(0.0001, unitReg.inch)
        iterated -= 1

    # After loop, set the closest diameter found
    OUT_FILM_C.Number(closest_diameter, CD_drill, Pressure_Chamber * Pressure_Drop_Fuel)
    # Initialize the actual values for the selected diameter
    if manual_override[0] == 0:
        Actual_hole_number = OUT_FILM_C.Hole_Number
    else:
        Actual_hole_number = manual_override[1]
    OUT_FILM_C.Actual(closest_diameter, Actual_hole_number)
    film_Cool_Doublet = [Original_Film_Cooling_Diameter, closest_diameter, drill_size, Actual_hole_number]
    return film_Cool_Doublet


def reinitialize_fuel(OUT_FILM_C, FUEL_CORE, OX_CORE, film_Cool_Doublet):
    Actual_percentage = OUT_FILM_C.Area_Actual / FUEL_CORE.Area_Actual
    OUT_FILM_C.mdot = (TotalMdot - OX_CORE.mdot) * Actual_percentage
    FUEL_CORE.mdot = TotalMdot - OX_CORE.mdot - OUT_FILM_C.mdot
    OUT_FILM_C.Actual(film_Cool_Doublet[1], film_Cool_Doublet[3])
    return OUT_FILM_C, FUEL_CORE, OX_CORE

def droplet_sizing(FUEL_CORE, OX_CORE, Doublet_Diameter_Fuel, Doublet_Diameter_LOX, viscosity_f, SurfaceTens_f):
    """Calculate droplet sizes using Dickerson and NASA methods."""

    # DROPLET SIZING
    # R. A. DlCKERSON method
    D_f_Dickerson = Q_( 1e5 * np.power(Doublet_Diameter_Fuel.magnitude, 0.27) *  np.power(Doublet_Diameter_LOX.magnitude, 0.023) / \
        (  np.power(FUEL_CORE.Velocity_Actual.magnitude, 0.74) * np.power(OX_CORE.Velocity_Actual.magnitude, 0.33)), unitReg.micron   )
    P_dynamic_OX = (0.5 * OX_CORE.rho * OX_CORE.Velocity_Actual**2).to(unitReg.psi)
    P_dynamic_FUEL = (0.5 * FUEL_CORE.rho * FUEL_CORE.Velocity_Actual**2).to(unitReg.psi)
    # Now use these in the NASA SP-8089 droplet size calculation
    P_D = (P_dynamic_OX / P_dynamic_FUEL).magnitude
    viscosity_shellwax = Q_(2.69e-3, unitReg.pound / (unitReg.foot * unitReg.second))  # lbm/(ft-sec)
    rho_shellwax = Q_(47.7, unitReg.pound / unitReg.foot**3)  # lbm/ft^3
    SurfaceTens_shellwax = Q_(17, unitReg.dyne / unitReg.centimeter)  # dynes/cm
    viscosity_f = viscosity_f.to(unitReg.pound / (unitReg.foot * unitReg.second))
    SurfaceTens_f = Q_(SurfaceTens_f.to(unitReg.dyne / unitReg.centimeter))
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
        unitReg.micron)
    return D_f_Dickerson, D_f_NASA


def calculate_vaporization_time_and_chamber_length(OX_CORE, FUEL_CORE, specific_heat_p_f, thermal_conductivity_f,
                                                   D_f_Dickerson, D_f_NASA):
    """Calculate vaporization time and chamber length for given droplet sizes."""
    B = ((Q_(5800, unitReg.degR) - Q_(525, unitReg.degF)) * Q_(0.475,  unitReg.BTU / unitReg.pound / unitReg.degR) / Q_(125,  unitReg.BTU / unitReg.pound )).to(unitReg.inch / unitReg.inch)
    #B = 26.8 # B for most hydrocarbon fuels are between 5 & 20.
    Vaporize_time_Dickerson = (FUEL_CORE.rho * specific_heat_p_f * (D_f_Dickerson/2)**2 / (2* thermal_conductivity_f * np.log(1 + B))).to(unitReg.second)
    Vaporize_time_NASA = (FUEL_CORE.rho * specific_heat_p_f * (D_f_NASA/2)**2 / (2* thermal_conductivity_f * np.log(1 + B.magnitude))).to(unitReg.second)
    gamma_lox_rad = np.deg2rad(OX_CORE.gamma.magnitude)
    gamma_fuel_rad = np.deg2rad(FUEL_CORE.gamma.magnitude)
    # Calculate momentum components
    momentum_x = (OX_CORE.mdot * OX_CORE.Velocity_Actual * np.cos(gamma_lox_rad) +
                  FUEL_CORE.mdot * FUEL_CORE.Velocity_Actual * np.cos(gamma_fuel_rad))
    momentum_y = (OX_CORE.mdot * OX_CORE.Velocity_Actual * np.sin(gamma_lox_rad) +
                  FUEL_CORE.mdot * FUEL_CORE.Velocity_Actual * np.sin(gamma_fuel_rad))
    Resultant_angle = np.arctan2(momentum_y,momentum_x)
    # Calculate total resultant impingement velocity from x and y components
    V_impingement = (np.sqrt(momentum_x**2 + momentum_y**2) / (OX_CORE.mdot + FUEL_CORE.mdot)).to(unitReg.feet / unitReg.second)
    ic(V_impingement)
    Travel_Length_Dickerson = (V_impingement *Vaporize_time_Dickerson).to(unitReg.inch)
    Travel_Length_NASA = (V_impingement *Vaporize_time_NASA).to(unitReg.inch)
    Chamber_Length_Dickerson = Travel_Length_Dickerson * np.cos(Resultant_angle)
    Chamber_Length_NASA = Travel_Length_NASA * np.cos(Resultant_angle)
    Dickerson = [D_f_Dickerson, Vaporize_time_Dickerson, Travel_Length_Dickerson, Chamber_Length_Dickerson]
    NASA = [D_f_NASA, Vaporize_time_NASA, Travel_Length_NASA, Chamber_Length_NASA]
    return Dickerson, NASA


def table_results(fuel_Doublet, film_Cool_Doublet, Dickerson, NASA, Doublet_Diameter_LOX, OX_CORE, FUEL_CORE,OUT_FILM_C,FuelName,gamma_FUEL_original):
    Original_Fuel_Diameter, closest_Fuel_Diameter, fuel_Doublet_Drill = fuel_Doublet
    Original_Film_Cooling_Diameter, closest_Film_Cool_Diameter, film_Cool_Doublet_Drill, _ = film_Cool_Doublet
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
    rounded_Table.add_row(["Fuel Diameter", f"{Original_Fuel_Diameter:.5f~}", closest_Fuel_Diameter],)
    rounded_Table.add_row(["Film Cooling ", f"{Original_Film_Cooling_Diameter:.5f~}", closest_Film_Cool_Diameter],)
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


def plot_results(OX_CORE, FUEL_CORE,OUT_FILM_C,Spacing,Rgamma_lox, ri, FilmCoolingSpacing, AimX, AimY, x_profile, y_profile, Chamber_ContourX, Chamber_ContourY):
    Points = 1000 #also used in other plots if this section ends up getting deleted

    # -------------- Constants and parameters -------------- #
    gamma_lox = OX_CORE.gamma.magnitude  # degrees and making the variables work below since i made the matlab version of this first and converted to python with ChatGPT
    gamma_fuel = FUEL_CORE.gamma.magnitude  # degrees


    # -------------- Plotting the aerospike nozzle contour -------------- #
    plt.figure()
    plt.plot(x_profile, y_profile, linewidth=2) 
    plt.axhline(0, linewidth=2, color="b")


    # -------------- Calculations and Plotting for impingement point and Resultant Point -------------- #
    x = Spacing.magnitude * np.sin(np.radians(90 + gamma_fuel)) / np.sin(np.radians(gamma_lox - gamma_fuel)) * np.sin(np.radians(90 - gamma_lox)) 
    y = Spacing.magnitude * np.sin(np.radians(90 + gamma_fuel)) / np.sin(np.radians(gamma_lox - gamma_fuel)) * np.cos(np.radians(90 - gamma_lox)) + Rgamma_lox.magnitude
    plt.plot(x, y, "o")
    yprime = AimY
    xprime = AimX
    plt.plot(xprime, yprime, "o")
    ic(x)
    ic(y)
    ic(xprime)
    ic(yprime)

    # -------------- Creating all linspaces needed for the following plots -------------- #
    x_graph = np.linspace(0, max(x_profile), Points)
    x_angled_lines = np.linspace(0, x, Points)  # Up to the impingement point for FUEL line
    #print(FilmCoolingSpacing[0]* np.tan( np.pi / 2  - np.radians(IN_FILM_C.gamma.magnitude)))
    x_filmcooling_outer = np.linspace(0, FilmCoolingSpacing.magnitude* np.tan( np.pi / 2  - np.abs(np.radians(OUT_FILM_C.gamma.magnitude))), Points)

    # -------------- Making and Plotting a Shitty Chamber Drawing  -------------- 
    plt.plot(Chamber_ContourX, Chamber_ContourY, "k", linewidth=2)

    # -------------- Plotting Horizontal Lines axial lines (References for angles) -------------- #
    plt.axhline(Rgamma_lox.magnitude, linestyle="--")
    plt.axhline((Rgamma_lox + Spacing).magnitude, linestyle="--")
    plt.axhline(y, linestyle="--") #need y before this plot. THis is the horizontal line to represent the axial POV from the resultant angle
    plt.axhline((ri - FilmCoolingSpacing).magnitude , linestyle="dotted")

    # -------------- Plotting the angled Prop Lines -------------- #
    gamma_lox_line = np.tan(np.radians(gamma_lox)) * x_angled_lines + Rgamma_lox.magnitude
    plt.plot(x_angled_lines, gamma_lox_line,"g")
    gamma_fuel_line = np.tan(np.radians(gamma_fuel)) * x_angled_lines + (Rgamma_lox + Spacing).magnitude
    plt.plot(x_angled_lines, gamma_fuel_line,"r")

    # -------------- Plotting the film cooling lines -------------- #
    outercooling_line = np.tan(np.radians(OUT_FILM_C.gamma.magnitude)) * x_filmcooling_outer + (ri - FilmCoolingSpacing).magnitude
    plt.plot(x_filmcooling_outer, outercooling_line,"r")


    # -------------- Plotting the resultant line -------------- #
    tan_resultant =     (OX_CORE.mdot.magnitude * OX_CORE.Velocity_Actual.magnitude * np.sin(np.deg2rad(OX_CORE.gamma.magnitude)) \
        + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_Actual.magnitude * np.sin(np.deg2rad(FUEL_CORE.gamma.magnitude)))/ \
        (OX_CORE.mdot.magnitude * OX_CORE.Velocity_Actual.magnitude * np.cos(np.deg2rad(OX_CORE.gamma.magnitude)) \
        + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity_Actual.magnitude * np.cos(np.deg2rad(FUEL_CORE.gamma.magnitude))) 
    resultant_y_intercept = y - tan_resultant * x
    ResultantX = x_graph / x_graph[-1] * (xprime - x) + x
    ResultantY = tan_resultant * ResultantX + resultant_y_intercept
    plt.plot(ResultantX, ResultantY, "y", linewidth=2)


    # -------------- Plotting the resultant error line And Error -------------- #
    a,b,c = -tan_resultant,1,-resultant_y_intercept
    Error = Q_(np.abs(a*xprime + b*yprime +c)/np.sqrt(a**2 + b**2),unitReg.inch)
    print(f'The errors from rounding to nearest drill size, and half degree for angles: resulted in missing the target by {Error :.4f~}')


    # -------------- Plotting the arcs that show the angles -------------- #
    arc_lox = patches.Arc((0, Rgamma_lox.magnitude), x, x, 
                          angle=0, theta1=0, theta2=gamma_lox, color='green', label=f'Gamma_LOX: {gamma_lox} deg')
    plt.gca().add_patch(arc_lox)
    arc_fuel = patches.Arc((0, (Rgamma_lox + Spacing).magnitude), x, x, 
                           angle=0, theta1=gamma_fuel, theta2=0, color='red', label=f'Gamma_FUEL: {gamma_fuel} deg')
    plt.gca().add_patch(arc_fuel)
    delta = Q_((np.arctan(tan_resultant)) *180 / np.pi, unitReg.degrees)
    arc_delta = patches.Arc((x, y), 3*x, 3*x, 
                           angle=0, theta1=0, theta2=delta.magnitude, color='y', label=f'Resultant_Angle: {delta} deg')
    plt.gca().add_patch(arc_delta)
    arc_filmcooling_outer = patches.Arc((0, (ri - FilmCoolingSpacing).magnitude), 1.5*x, 1.5*x, 
                          angle=0, theta1=0, theta2=OUT_FILM_C.gamma.magnitude, color='red', label=f'Gamma_FilmCooling_Outer: {OUT_FILM_C.gamma.magnitude} deg')
    plt.gca().add_patch(arc_filmcooling_outer)


    # -------------- Plotting the angles onto the graph for easier readibility -------------- #
    x_offset = -x/2  # Adjust as needed
    y_offset_lox = -3*(y - Rgamma_lox.magnitude)/4  
    y_offset_fuel = 5*((Spacing + Rgamma_lox).magnitude - y)/8  
    x_offset_delta = 3*(xprime-x)/4  
    y_offset_delta = 5*(yprime-y)/8
    x_offset_Filmcooling = x
    y_offset_Filmcooling_Outer =  (FilmCoolingSpacing/4).magnitude
    FUEL_CORE_gamma_writing = abs(FUEL_CORE.gamma.magnitude)
    plt.text(x + x_offset,y + y_offset_lox, rf'$\gamma_{{\mathrm{{OX}}}}: {OX_CORE.gamma.magnitude:.2f}^\circ$', color='green')
    plt.text(x + x_offset,y + y_offset_fuel, rf'$\gamma_{{\mathrm{{F}}}}: {FUEL_CORE_gamma_writing:.2f}^\circ$', color='red')
    plt.text(x + x_offset_delta, y + y_offset_delta, rf'$\delta: {delta.magnitude:.2f}^\circ$', color='y')
    plt.text(x_offset_Filmcooling, (ri - FilmCoolingSpacing).magnitude + y_offset_Filmcooling_Outer, rf'$\gamma_{{\mathrm{{FCO}}}}: {OUT_FILM_C.gamma.magnitude:.2f}^\circ$', color='red')


    # -------------- Extra Plotting Shit -------------- #
    plt.legend(['Spike Contour', 'Centerline', 'Impingement Point', 'Aim Point','Chamber Contour', 'Gamma_(OX) Straight Line', 'Gamma_(FUEL) Straight Line',
                'Resultant Straight Line', 'Film Cooling Outer Straight Line',
                  f'Gamma_(OX) Angled Line {OX_CORE.gamma :.3f~}', 
                 f'Gamma_(FUEL) Angled Line {FUEL_CORE.gamma :.3f~}', 
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