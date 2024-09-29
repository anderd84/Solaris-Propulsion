from matplotlib import pyplot as plt, patches
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from pint import UnitRegistry
import numpy as np
import pyromat as pm


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



# -------------- Pyromat Shit -------------- #
pm.config['unit_pressure'] = 'psi'
pm.config['unit_mass'] = 'lbm'
pm.config['unit_matter'] = 'lbm'
pm.config['unit_length'] = 'ft'
pm.config['unit_volume'] = 'ft3'
N2 = pm.get('mp.N2')
LOX = pm.get('mp.O2')
CO2 = pm.get('mp.CO2')

# -------------- Unit Shit -------------- #
ureg = UnitRegistry()
ureg.default_system = 'US'
Q_ = ureg.Quantity



# -------------- Constants -------------- #
Prescott_pressure = Q_(12.04, ureg.force_pound / ureg.inch**2)
CD_drill = 0.7 #Constant for Sharp Edged Orifices         
g0 = Q_(32.174, ureg.foot / ureg.second**2)
JetADensity = 51.15666

# -------------- Lox Dewar Pressure -------------- #
def LOXDensity(Lox_Dewar_Pressure):
    Lox_Absolute_Pressure = Q_(Lox_Dewar_Pressure, ureg.force_pound / ureg.inch**2) + Prescott_pressure
    LOX_Sat_Temp = LOX.Ts(p=Lox_Absolute_Pressure.magnitude )
    LOX_Sat_Dens = LOX.ds(T=LOX_Sat_Temp)
    LOX_Sat_Dens = LOX_Sat_Dens[0][0]
    return LOX_Sat_Dens

# -------------- ROT -------------- #
Pressure_Drop_Fuel = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Drop_Lox = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)



# -------------- $ 4 Different PROP FLOWS -------------- #
def PROPFLOWS(mdots,Film_Cooling,gammas,Lox_Dewar_Pressure):
    OX_CORE = PROP(gamma=gammas[0], mdot=mdots[0], rho=LOXDensity(Lox_Dewar_Pressure))
    FUEL_CORE = PROP(gamma = gammas[1], mdot = mdots[1], rho=JetADensity) #gamma zero for this one because it's the initialized guess just making the FUEL CORE class requires it ( should change when moving to data classes)
    OUT_FILM_C = PROP(gamma = gammas[2], mdot = Film_Cooling[0]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
    IN_FILM_C = PROP(gamma = gammas[3], mdot = Film_Cooling[1]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
    return OX_CORE,FUEL_CORE,OUT_FILM_C,IN_FILM_C
