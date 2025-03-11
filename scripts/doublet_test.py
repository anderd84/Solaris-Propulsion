import numpy as np
from icecream import ic

from fluids.fluid import CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox
import general.design as DESIGN
from general.units import Q_, unitReg
from injector import doublet
from nozzle import plug

# davids stuff
Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)

plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))

plugx = np.array([p.x for p in plugC])
plugr = np.array([p.r for p in plugC])
offset = min(plugx)
cowlx = np.array([p.x for p in cowlC])
cowlr = np.array([p.r for p in cowlC])


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
AimX = straightLength.m_as("in")
AimY = aimpoint[1]
chamberatInjectorRadius = DESIGN.chamberatInjectorRadius
OX_CORE,FUEL_CORE,OUT_FILM_C,viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f = doublet.initialize_prop_flows(
    Film_Cooling, oxImpingeAngle, fuelInitalImpingeAngle, filmImpingeAngle,Lox_Dewar_Pressure, CD_drill, Pressure_Drop_Lox,  Pressure_Drop_Fuel, AirTemperature, AirPressure, FuelName, Pressure_Chamber)
gamma_FUEL_original, FUEL_CORE = doublet.calculate_fuel_impinge_angle(OX_CORE, FUEL_CORE, AimY, Rgamma_lox, Spacing, AimX)
fuel_Doublet = doublet.calculate_fuel_diameters(OX_CORE, FUEL_CORE, Pressure_Chamber, Doublet_Diameter_LOX, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox)
manual_override = [0 , 51] # wheteher you actual override the cell or not, The new hole number if you need to
film_Cool_Doublet= doublet.calculate_film_cooling_diameters(OUT_FILM_C, FilmCoolingSpacing, Pressure_Chamber, CD_drill, Pressure_Drop_Fuel, chamberatInjectorRadius,manual_override)
ic(OUT_FILM_C.Area_Actual / FUEL_CORE.Area_Actual)


OUT_FILM_C, FUEL_CORE, OX_CORE = doublet.reinitialize_fuel(OUT_FILM_C, FUEL_CORE, OX_CORE, film_Cool_Doublet)

#fuel_Doublet = doublet.calculate_fuel_diameters(OX_CORE, FUEL_CORE, Pressure_Chamber, Doublet_Diameter_LOX, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox)
D_f_Dickerson, D_f_NASA =doublet.droplet_sizing(FUEL_CORE, OX_CORE, fuel_Doublet[1], Doublet_Diameter_LOX, viscosity_f, SurfaceTens_f)
Dickerson, NASA = doublet.calculate_vaporization_time_and_chamber_length(OX_CORE, FUEL_CORE, specific_heat_p_f, thermal_conductivity_f,D_f_Dickerson, D_f_NASA)
doublet.table_results(fuel_Doublet, film_Cool_Doublet, Dickerson, NASA, Doublet_Diameter_LOX, OX_CORE, FUEL_CORE,OUT_FILM_C,FuelName,gamma_FUEL_original)

doublet.plot_results(OX_CORE, FUEL_CORE,OUT_FILM_C, Spacing, Rgamma_lox, ri, FilmCoolingSpacing, AimX, AimY, plugx - offset, plugr, cowlx - offset, cowlr)