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


# -------------- Design Inputs and Constants -------------- #    
#Constants
CD_drill = 0.7 #Constant for Sharp Edged Orifices         
g0 = Q_(32.174, ureg.foot / ureg.second**2)
Prescott_pressure = Q_(12.04, ureg.force_pound / ureg.inch**2)
#Design Input
mdots = np.array([5.29, 2.21])/4 #LOX_CORE, FUEL_CORE
Film_Cooling = np.array([0.05, 0.05]) #Outer Film Cooling Percentage, Inner Film Cooling Percentage
di = 6.5 #Internal Diameter of Chamber
ri = di / 2 #Internal Radius of Chamber
Spacing = 0.55  #Spacing between center of impingement Holes
Rgamma_lox = 1.65  #Radial distance between centerline and LOX hole
FilmCoolingSpacing = np.array([.25, .25]) #inches Inner, Outer
Pressure_Drop_Fuel = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Drop_Lox = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Chamber = Q_(300, ureg.force_pound / ureg.inch**2) #Chamber Pressure Pretty Obvious
Doublet_Diameter_LOX = Q_(0.125, ureg.inch)  #Design choise for DOublet Diameter size (Need to look more properly into it as 1/4 holes might make vaporization time too long)\
Lox_Dewar_Pressure = Q_(22, ureg.force_pound / ureg.inch**2) + Prescott_pressure

# -------------- Prop Initialization -------------- #
LOX_Sat_Temp = LOX.Ts(p=Lox_Dewar_Pressure.magnitude )
LOX_Sat_Dens = LOX.ds(T=LOX_Sat_Temp)
OX_CORE = PROP(gamma=25., mdot=mdots[0], rho=LOX_Sat_Dens[0][0])
FUEL_CORE = PROP(gamma = 0, mdot = mdots[1], rho=51.15666) #gamma zero for this one because it's the initialized guess just making the FUEL CORE class requires it ( should change when moving to data classes)
OUT_FILM_C = PROP(gamma = 30, mdot = Film_Cooling[0]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
IN_FILM_C = PROP(gamma = -30, mdot = Film_Cooling[1]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)