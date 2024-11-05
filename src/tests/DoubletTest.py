from matplotlib import pyplot as plt, patches
from Injector.Doublet_Functions import spike_contour
from Injector.InjectorCad import injector_cad_write
from Injector.Drill import drill_approximation
from Injector import ActualDoublet as doublet
from scipy.optimize import fsolve
from fluids.fluid import CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox
import numpy as np
from icecream import ic
import General.design as DESIGN


# -------------- Design Inputs -------------- #    
ri = DESIGN.chamberInternalRadius #Internal Diameter of Chamber
Spacing = DESIGN.Spacing  #Spacing between centear of impingement Holes
Rgamma_lox = DESIGN.oxHoleRadius  #Radial distance between centerline and LOX hole
Film_Cooling = DESIGN.percentFilmCooling #Outer Film Cooling Percentage
FilmCoolingSpacing = DESIGN.filmCoolingSpacing #inches Outer
Pressure_Chamber = DESIGN.chamberPressure #Chamber Pressure Pretty Obvious
Doublet_Diameter_LOX = DESIGN.oxDoubletDiameter  #Design choise for DOublet Diameter size (Need to look more properly into it as 1/4 holes might make vaporization time too long)\
oxImpingeAngle, fuelInitalImpingeAngle, filmImpingeAngle = DESIGN.oxImpingeAngle,0,DESIGN.filmImpingeAngle #Lox, fuel, outer FC  #All Angles from axial
Lox_Dewar_Pressure = DESIGN.oxDewarPressure
AirTemperature = DESIGN.prescottAmbientTemp
AirPressure = DESIGN.prescottAmbientPressure
FuelName = DESIGN.fuelName
Points = 1000 #also used in other plots if this section ends up getting deleted
AimX = 8
AimY = 3.66

OX_CORE,FUEL_CORE,OUT_FILM_C,viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f = doublet.initialize_prop_flows(
    Film_Cooling, oxImpingeAngle, fuelInitalImpingeAngle, filmImpingeAngle,Lox_Dewar_Pressure, CD_drill, Pressure_Drop_Lox,  Pressure_Drop_Fuel, AirTemperature, AirPressure, FuelName, Pressure_Chamber)
gamma_FUEL_original, FUEL_CORE = doublet.calculate_fuel_impinge_angle(OX_CORE, FUEL_CORE, AimY, Rgamma_lox, Spacing, AimX)
fuel_Doublet = doublet.calculate_fuel_diameters(OX_CORE, FUEL_CORE, Pressure_Chamber, Doublet_Diameter_LOX, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox)
film_Cool_Doublet = doublet.calculate_film_cooling_diameters(OUT_FILM_C, FilmCoolingSpacing, Pressure_Chamber, CD_drill, Pressure_Drop_Fuel, ri)
D_f_Dickerson, D_f_NASA =doublet.droplet_sizing(FUEL_CORE, OX_CORE, fuel_Doublet[1], Doublet_Diameter_LOX, viscosity_f, SurfaceTens_f)
Dickerson, NASA = doublet.calculate_vaporization_time_and_chamber_length(OX_CORE, FUEL_CORE, specific_heat_p_f, thermal_conductivity_f,D_f_Dickerson, D_f_NASA)
doublet.table_results(fuel_Doublet, film_Cool_Doublet, Dickerson, NASA, Doublet_Diameter_LOX, OX_CORE, FUEL_CORE,OUT_FILM_C,FuelName,gamma_FUEL_original)

x_profile,y_profile = spike_contour(Points) #Will be the spike Contour

x_graph = np.linspace(0, max(x_profile), Points) # -------------- Making and Plotting a Shitty Chamber Drawing  -------------- 
thetaRange4 = np.linspace(np.radians(90), 0, Points)
Chamber_Cowl_r = 0.5  # in
Past_Peak = 1.15 #some terrible constant for shitty chamber I made drawn
ChamberX = x_graph / x_graph[-1] * AimX * Past_Peak
ChamberY = np.ones(len(ChamberX)) * ri.magnitude
ChamberArcX = ChamberX[-1] + Chamber_Cowl_r * np.cos(thetaRange4)
ChamberArcY = ChamberY[-1] + Chamber_Cowl_r * (-1 + np.sin(thetaRange4))
Chamber_ContourX = np.concatenate([ChamberX, ChamberArcX])
Chamber_ContourY = np.concatenate([ChamberY, ChamberArcY])

doublet.plot_results(OX_CORE, FUEL_CORE,OUT_FILM_C,Spacing,Rgamma_lox, ri, FilmCoolingSpacing, AimX, AimY, x_profile, y_profile, Chamber_ContourX, Chamber_ContourY)