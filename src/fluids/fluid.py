from matplotlib import pyplot as plt, patches
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from pint import UnitRegistry
import numpy as np
import pyromat as pm
import CoolProp.CoolProp as CP

OF = 2.4
TotalMdot = 9


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
def PROPFLOWS(Film_Cooling,gammas,Lox_Dewar_Pressure):
    mdots = np.array([OF*TotalMdot/(1+OF), TotalMdot/(1+OF)]) #LOX_CORE, FUEL_CORE lbm/s
    OX_CORE = PROP(gamma=gammas[0], mdot=mdots[0], rho=LOXDensity(Lox_Dewar_Pressure))
    FUEL_CORE = PROP(gamma = gammas[1], mdot = mdots[1]*(1 - Film_Cooling[0] - Film_Cooling[1]), rho=JetADensity) #gamma zero for this one because it's the initialized guess just making the FUEL CORE class requires it ( should change when moving to data classes)
    OUT_FILM_C = PROP(gamma = gammas[2], mdot = Film_Cooling[0]* mdots[1], rho = FUEL_CORE.rho)
    IN_FILM_C = PROP(gamma = gammas[3], mdot = Film_Cooling[1]* mdots[1], rho = FUEL_CORE.rho)
    print("checking that mdots were calculated right... error =", 
          Q_(mdots[0], ureg.pound / ureg.second) + Q_(mdots[1], ureg.pound / ureg.second) - OX_CORE.mdot - FUEL_CORE.mdot - OUT_FILM_C.mdot - IN_FILM_C.mdot)
        


    return OX_CORE,FUEL_CORE,OUT_FILM_C,IN_FILM_C

# -------------- $ Cooling Equations -------------- #
#def bartzConv()


# -------------- $ Kerosene Properties -------------- #
def get_fluid_properties(name, temperature_R, pressure_psi):
    temperature = Q_(temperature_R, ureg.degR)
    pressure = Q_(pressure_psi, ureg.psi)
    temperature_SI = temperature.to(ureg.kelvin).magnitude
    pressure_SI = pressure.to(ureg.pascal).magnitude
    
    # Viscosity
    viscosity_SI = Q_(CP.PropsSI('V', 'T', temperature_SI, 'P', pressure_SI, name), ureg.pascal * ureg.second)
    viscosity = viscosity_SI.to(ureg.pound / ureg.foot / ureg.second)
    
    # Specific heat at constant pressure
    specific_heat_p_SI = Q_(CP.PropsSI('C', 'T', temperature_SI, 'P', pressure_SI, name), ureg.joule / ureg.kilogram / ureg.kelvin)
    specific_heat_p = specific_heat_p_SI.to(ureg.BTU / ureg.pound / ureg.degR)
    
    # Specific heat at constant volume
    specific_heat_v_SI = Q_(CP.PropsSI('O', 'T', temperature_SI, 'P', pressure_SI, name), ureg.joule / ureg.kilogram / ureg.kelvin)
    
    # Gamma (Cp/Cv)
    gamma = specific_heat_p_SI / specific_heat_v_SI
    
    # Thermal conductivity
    thermal_conductivity_SI = Q_(CP.PropsSI('L', 'T', temperature_SI, 'P', pressure_SI, name), ureg.watt / ureg.meter / ureg.kelvin)
    thermal_conductivity = thermal_conductivity_SI.to(ureg.BTU / ureg.foot / ureg.hour / ureg.degR)
    
    # Density (kg/m³)
    density_SI = Q_(CP.PropsSI('D', 'T', temperature_SI, 'P', pressure_SI, name), ureg.kilogram / ureg.meter**3)
    density = density_SI.to(ureg.pound / ureg.foot**3)  # Convert to lbm/ft³
    
    # Prandtl number (dimensionless)
    prandtl = CP.PropsSI('PRANDTL', 'T', temperature_SI, 'P', pressure_SI, name)
    
    # Quality (dimensionless) - vapor mass fraction (0 = saturated liquid, 1 = saturated vapor)
    quality = CP.PropsSI('Q', 'T', temperature_SI, 'P', pressure_SI, name)
    
    # Phase (integer: 0 = unknown, 1 = liquid, 2 = vapor, 3 = supercritical)
    phase = CP.PropsSI('PHASE', 'T', temperature_SI, 'P', pressure_SI, name)
    
    # Coefficient of thermal expansion (1/K)
    alpha_SI = Q_(CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', temperature_SI, 'P', pressure_SI, name), 1 / ureg.kelvin)
    
    # Thermal diffusivity (m²/s)
    thermal_diffusivity_SI = (thermal_conductivity_SI / (density_SI * specific_heat_p_SI)).to(ureg.meter**2 / ureg.second)
    thermal_diffusivity = thermal_diffusivity_SI.to(ureg.foot**2 / ureg.second)  # Convert to ft²/s
    
    # Return all the properties
    return (viscosity, specific_heat_p, gamma, thermal_conductivity, density, prandtl, 
            quality, phase, alpha_SI, thermal_diffusivity)
