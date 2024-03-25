import numpy as np
import scipy as sp
from scipy import linalg
rg = np.random.default_rng()
import matplotlib.pyplot as plt
from scipy import integrate
from pint import UnitRegistry
ureg = UnitRegistry()
ureg.default_system = 'US'
Q_ = ureg.Quantity
class PROP:
    def __init__(self, gamma, mdot, rho):
        self.gamma = Q_(gamma, ureg.degrees)
        self.mdot = Q_(mdot, ureg.pound / ureg.second)
        self.rho = Q_(rho, ureg.pound / ureg.foot**3)
    def Velocity(self,Cd , Pressure_Diff):
        Velocity = Cd * np.sqrt(2 * Pressure_Diff/ self.rho *gc)
        return Velocity.to(ureg.feet / ureg.second)
    def Area(self, Cd, Pressure_Diff):
        Area = (self.mdot / (Cd * np.sqrt(2*self.rho* Pressure_Diff*gc)))
        return Area.to(ureg.inch**2)
    
#Constants
CD_drill = 0.7 #Constant for Sharp Edged Orifices
gc = Q_(32.174,  ureg.pound * ureg.ft / ( ureg.force_pound * ureg.second**2)) 
#(lb·ft)/(lbf·s2) #might not be needed Apparent Pint library takes care of it

#Design Input
#gamma defined as angle after faceplate off the axial direction. Positive defined as away from center of Aerospike
mdots = np.array([5.29, 2.21]) #LOX_CORE, FUEL_CORE
Film_Cooling = np.array([0.08, 0.08]) #Outer Film Cooling, Inner Film Cooling
Pressure_Drop_Fuel = 0.2
Pressure_Drop_Lox = 0.2
Pressure_Chamber = Q_(300, ureg.force_pound / ureg.inch**2)
OX_CORE = PROP(gamma=-15, mdot=mdots[0], rho=56.794)
FUEL_CORE = PROP(gamma = 12, mdot = mdots[1], rho=51.15666)
OUT_FILM_C = PROP(gamma = 30, mdot = Film_Cooling[0]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
IN_FILM_C = PROP(gamma = -30, mdot = Film_Cooling[1]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)


print(f"Total OX Doublets Velocity: {OX_CORE.Velocity(CD_drill, Pressure_Chamber * (Pressure_Drop_Lox)):.3f~}")
print(f"Total OX Orifice Area Doublets: {OX_CORE.Area(CD_drill, Pressure_Chamber * (Pressure_Drop_Lox)):.3f~}")
FUEL_mdot_tot = FUEL_CORE.mdot + OUT_FILM_C.mdot + IN_FILM_C.mdot
print(f"For Total Fuel mdot in normal units: {FUEL_mdot_tot:.3f~}")
print(f"For Total Fuel mdot in commie units: {FUEL_mdot_tot.to(ureg.kilogram / ureg.second):.3f~}")
