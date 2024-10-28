from matplotlib import pyplot as plt, patches
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from pint import UnitRegistry
import numpy as np
import pyromat as pm
import CoolProp.CoolProp as CP
from icecream import ic
from rocketprops.rocket_prop import get_prop

# -------------- Pyromat Shit -------------- #
pm.config['unit_pressure'] = 'psi'
pm.config['unit_mass'] = 'lbm'
pm.config['unit_matter'] = 'lbm'
pm.config['unit_length'] = 'ft'
pm.config['unit_volume'] = 'ft3'
pm.config['unit_volume'] = 'ft3'
pm.config['unit_temperature'] = 'Rankine'
N2 = pm.get('mp.N2')
LOX = pm.get('mp.O2')
CO2 = pm.get('mp.CO2')

# -------------- Unit Shit -------------- #
ureg = UnitRegistry()
ureg.default_system = 'US'
Q_ = ureg.Quantity
ureg.default_format = "~P"  # Compact unit formatting
ureg.define("lbmol = 453.59237 * mol")  # 1 lbmol = 453.59237 mol (since 1 lb = 453.59237 g)

OF = 2.0
TotalMdot = 8
Chamber_Press = 300.0 #psia
from rocketcea.cea_obj import CEA_Obj; Combustion=CEA_Obj(oxName='LOX', fuelName='RP-1');
Combustion_Temp = Q_(Combustion.get_Tcomb(Pc=Chamber_Press, MR=OF), ureg.degR)
IspVac, Cstar, Tc, MW, gamma = Combustion.get_IvacCstrTc_ThtMwGam(Pc=Chamber_Press, MR=OF, eps=7.11)
MW = Q_(MW, ureg.pound / ureg.lbmol)
R_universal = Q_(10.731577089016, ureg.psi * ureg.foot**3 / (ureg.lbmol * ureg.degR))  # Universal gas constant in ft·lbf/(lbmol·°R)
R_specific = (R_universal / MW).to(ureg.foot * ureg.pound_force / (ureg.pound * ureg.degR))
#ic(R_specific)


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
        self.Number_Actual = Number
        return Velocity_Actual, Area_Actual           





# -------------- Constants -------------- #
Prescott_pressure = Q_(12.04, ureg.force_pound / ureg.inch**2)
CD_drill = 0.7 #Constant for Sharp Edged Orifices         
g0 = Q_(32.174, ureg.foot / ureg.second**2)
JetADensity = 51.15666

# -------------- Lox Dewar Pressure -------------- #
def LOXDensity(Lox_Dewar_Pressure):
    LOX_Absolute_Pressure = Q_(Lox_Dewar_Pressure, ureg.force_pound / ureg.inch**2) + Prescott_pressure
    LOX_Sat_Temp = LOX.Ts(p=LOX_Absolute_Pressure.magnitude )
    LOX_Sat_Dens = LOX.ds(T=LOX_Sat_Temp)
    LOX_Sat_Dens = LOX_Sat_Dens[0][0]
    return LOX_Sat_Dens

# -------------- ROT -------------- #
Pressure_Drop_Fuel = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Drop_Lox = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)



# -------------- $ 4 Different PROP FLOWS -------------- #
def PROPFLOWS(Film_Cooling,gammas,Lox_Dewar_Pressure, AirTemperature, AirPressure,fuel_name):
    temperature_R = AirTemperature.to(ureg.degR)
    pressure_psi = AirPressure.to(ureg.psi)


    properties_fuel = get_fluid_properties(fuel_name,temperature_R.magnitude, pressure_psi.magnitude)

    (viscosity_f, specific_heat_p_f, gamma_f, thermal_conductivity_f, density_f, prandtl_f, alpha_f, thermal_diffusivity_f, SurfaceTens_f) = properties_fuel

    mdots = np.array([OF*TotalMdot/(1+OF), TotalMdot/(1+OF)]) #LOX_CORE, FUEL_CORE lbm/s
    OX_CORE = PROP(gamma=gammas[0], mdot=mdots[0], rho=LOXDensity(Lox_Dewar_Pressure))
    FUEL_CORE = PROP(gamma = gammas[1], mdot = mdots[1]*(1 - Film_Cooling[0] - Film_Cooling[1]), rho=density_f) #gamma zero for this one because it's the initialized guess just making the FUEL CORE class requires it ( should change when moving to data classes)
    OUT_FILM_C = PROP(gamma = gammas[2], mdot = Film_Cooling[0]* mdots[1], rho = FUEL_CORE.rho)
    IN_FILM_C = PROP(gamma = gammas[3], mdot = Film_Cooling[1]* mdots[1], rho = FUEL_CORE.rho)
 
#    print("checking that mdots were calculated right... error =", 
#          Q_(mdots[0], ureg.pound / ureg.second) + Q_(mdots[1], ureg.pound / ureg.second) - OX_CORE.mdot - FUEL_CORE.mdot - OUT_FILM_C.mdot - IN_FILM_C.mdot)
        


    return OX_CORE,FUEL_CORE,OUT_FILM_C,IN_FILM_C, viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f





# -------------- $ Cooling Equations -------------- #
#def bartzConv()
#rp1 = get_prop("RP1")
#rp1.summ_print()


# Function to retrieve fluid properties using RocketProps defaults (English Engineering units)
def get_fluid_properties(name, temperature_R, pressure_psi):
    # Ensure inputs are in the correct units for RocketProps methods
    temperature = Q_(temperature_R, ureg.degR)
    pressure = Q_(pressure_psi, ureg.psi)
    
    # Create an fluid object
    fluid = get_prop(name)
    
    # Retrieve dynamic properties from RocketProps at the given temperature and pressure
    viscosity_poise = Q_(fluid.ViscAtTdegR(temperature.magnitude), ureg.poise)  # Viscosity in poise
    viscosity = viscosity_poise.to(ureg.pound / ureg.foot / ureg.second)  # Convert to lb/ft·s using .to()
    
    specific_heat_p = Q_(fluid.CpAtTdegR(temperature.magnitude), ureg.BTU / ureg.pound / ureg.degR)  # Cp in BTU/lbm-R
    thermal_conductivity = Q_(fluid.CondAtTdegR(temperature.magnitude), ureg.BTU / ureg.foot / ureg.hour / ureg.degR)  # Thermal conductivity
    
    # Density (lbm/ft³) using compressed specific gravity and accurate water density
    specific_gravity = fluid.SG_compressed(temperature.magnitude, pressure.magnitude)  # Specific gravity in g/ml
    water_density = Q_(1, ureg.gram / ureg.milliliter).to(ureg.pound / ureg.foot**3)  # Density of water in lbm/ft³
    density = specific_gravity * water_density  # Convert specific gravity to density for RP-1 in lbm/ft³
    
    # Calculate Prandtl number (dimensionless)
    prandtl = (viscosity * specific_heat_p / thermal_conductivity).to('dimensionless')
    
    # Estimate gamma (Cp/Cv) assuming Cv ~ Cp/1.25 if no direct Cv is available
    specific_heat_v = specific_heat_p / 1.25
    gamma = specific_heat_p / specific_heat_v
    
    # Coefficient of thermal expansion (1/°R) as an approximation
    alpha = Q_(1 / (temperature.magnitude * 60), 1 / ureg.degR)
    
    # Thermal diffusivity (ft²/s)
    thermal_diffusivity = thermal_conductivity / (density * specific_heat_p) * ureg.foot**2 / ureg.second

    # Surface Tension (lbf/ft)
    SurfaceTens = Q_(fluid.SurfAtTdegR(temperature.magnitude), ureg.pound_force / ureg.foot)
    
    # Return all properties in English Engineering units
    return (viscosity, specific_heat_p, gamma, thermal_conductivity, density, prandtl, alpha, thermal_diffusivity, SurfaceTens)

