import numpy as np
import pyromat as pm
import rocketprops 
from rocketprops.rocket_prop import get_prop
from scipy.optimize import fsolve

from general.units import Q_, unitReg
import general.design as DESIGN

# -------------- Pyromat Shit -------------- #
pm.config['unit_pressure'] = 'psi'
pm.config['unit_mass'] = 'lbm'
pm.config['unit_matter'] = 'lbm'
pm.config['unit_length'] = 'ft'
pm.config['unit_volume'] = 'ft3'
pm.config['unit_volume'] = 'ft3'
pm.config['unit_temperature'] = 'Rankine'
LOX = pm.get('mp.O2')

OF = DESIGN.OFratio
TotalMdot = DESIGN.totalmdot
OxMdot = DESIGN.Oxidizer_Total
FuelMdot = DESIGN.Fuel_Total
Chamber_Press = DESIGN.chamberPressure


class PROP:
    def __init__(self, gamma, mdot, rho):
        """This is teh __init__ fn that will do a thing
        Args:
            gamma (float): Angle off axial direction with positive away from centerbody
            mdot (float): mdot of prop
            rho (float): density of prop
        """        
        self.gamma = gamma
        self.mdot = mdot
        self.rho = rho


    def Velocity(self,Cd , Pressure_Diff) -> float:
        """Velocity Function using Dr. Whites Eq 16.5
        Args:
            Cd (float): Coefficient of Discharge
            Pressure_Diff (float): Pressure difference through injector (Estimated metric) in psi
        Returns:
            float: Velocity Out of the injector parallel to flow
        """
        Velocity_init = (Cd * np.sqrt(2 * Pressure_Diff/ self.rho )).to(unitReg.feet / unitReg.second)
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
        Area_init = (self.mdot / (Cd * np.sqrt(2*self.rho* Pressure_Diff ))).to(unitReg.inch**2)
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
        Hole_Area = (np.pi * Hole_Diameter**2 /4).to(unitReg.inch**2)
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
        Area_Actual = (Number * 0.25 * Hole_Diameter**2 * np.pi).to(unitReg.inch**2)
        Velocity_Actual = (self.mdot/(self.rho * Area_Actual)).to(unitReg.feet / unitReg.second)
        
        self.Area_Actual = Area_Actual
        self.Velocity_Actual = Velocity_Actual
        self.Number_Actual = Number
        return Velocity_Actual, Area_Actual           





# -------------- Constants -------------- #
Prescott_pressure = DESIGN.prescottAmbientPressure
CD_drill = 0.7 #Constant for Sharp Edged Orifices         
g0 = Q_(32.174, unitReg.foot / unitReg.second**2)


# -------------- Lox Dewar Pressure -------------- #
def LOXDensity(Lox_Dewar_Pressure):
    LOX_Absolute_Pressure = Lox_Dewar_Pressure + Prescott_pressure
    LOX_Sat_Temp = LOX.Ts(p=LOX_Absolute_Pressure.magnitude )
    LOX_Sat_Dens = LOX.ds(T=LOX_Sat_Temp)
    LOX_Sat_Dens = LOX_Sat_Dens[0][0]
    return LOX_Sat_Dens

# -------------- ROT -------------- #
Pressure_Drop_Fuel = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)
Pressure_Drop_Lox = 0.2 #Pressure drop Percentage (ROT: Always in terms of Chamber Pressure)



# -------------- $ 4 Different PROP FLOWS -------------- #
def PROPFLOWS(Film_Cooling,oxImpingeAngle, fuelInitalImpingeAngle, filmImpingeAngle,Lox_Dewar_Pressure, AirTemperature, AirPressure,fuel_name):
    temperature_R = AirTemperature.to(unitReg.degR)
    pressure_psi = AirPressure.to(unitReg.psi)
    properties_fuel = get_fluid_properties(fuel_name,temperature_R, pressure_psi)
    (viscosity_f, specific_heat_p_f, gamma_f, thermal_conductivity_f, density_f, prandtl_f, alpha_f, thermal_diffusivity_f, SurfaceTens_f) = properties_fuel
    OX_CORE = PROP(gamma=oxImpingeAngle, mdot=OxMdot, rho=Q_(LOXDensity(Lox_Dewar_Pressure), unitReg.pound / unitReg.foot**3))
    FUEL_CORE = PROP(gamma = fuelInitalImpingeAngle, mdot = FuelMdot*(1 - Film_Cooling), rho=density_f) #gamma zero for this one because it's the initialized guess just making the FUEL CORE class requires it ( should change when moving to data classes)
    OUT_FILM_C = PROP(gamma =  filmImpingeAngle, mdot = Film_Cooling* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
 
#    print("checking that mdots were calculated right... error =", 
#          Q_(mdots[0], unitReg.pound / unitReg.second) + Q_(mdots[1], unitReg.pound / unitReg.second) - OX_CORE.mdot - FUEL_CORE.mdot - OUT_FILM_C.mdot - IN_FILM_C.mdot)
        


    return OX_CORE,FUEL_CORE,OUT_FILM_C,viscosity_f, specific_heat_p_f, thermal_conductivity_f, SurfaceTens_f

def get_cp(name, temperature_R):
    temperature = Q_(temperature_R.magnitude, unitReg.degR)
    fluid = get_prop(name)
    specific_heat_p = Q_(fluid.CpAtTdegR(temperature.magnitude), unitReg.BTU / unitReg.pound / unitReg.degR)
    return specific_heat_p

# Function to retrieve fluid properties using RocketProps defaults (English Engineering units)
def get_fluid_properties(name, temperature_R, pressure_psi):
    # Ensure inputs are in the correct units for RocketProps methods
    temperature = Q_(temperature_R.magnitude, unitReg.degR)
    
    pressure = Q_(pressure_psi.magnitude, unitReg.psi)
    
    # Create an fluid object
    fluid = get_prop(name)
    
    # Retrieve dynamic properties from RocketProps at the given temperature and pressure
    viscosity_poise = Q_(fluid.ViscAtTdegR(temperature.magnitude), unitReg.poise)  # Viscosity in poise
    viscosity = viscosity_poise.to(unitReg.pound / unitReg.foot / unitReg.second)  # Convert to lb/ft·s using .to()
    
    specific_heat_p = Q_(fluid.CpAtTdegR(temperature.magnitude), unitReg.BTU / unitReg.pound / unitReg.degR)  # Cp in BTU/lbm-R
    thermal_conductivity = Q_(fluid.CondAtTdegR(temperature.magnitude), unitReg.BTU / unitReg.foot / unitReg.hour / unitReg.degR)  # Thermal conductivity
    
    # Density (lbm/ft³) using compressed specific gravity and accurate water density
    specific_gravity = fluid.SG_compressed(temperature.magnitude, pressure.magnitude)  # Specific gravity in g/ml
    water_density = Q_(1, unitReg.gram / unitReg.milliliter).to(unitReg.pound / unitReg.foot**3)  # Density of water in lbm/ft³
    density = specific_gravity * water_density  # Convert specific gravity to density for RP-1 in lbm/ft³
    
    # Calculate Prandtl number (dimensionless)
    prandtl = (viscosity * specific_heat_p / thermal_conductivity).to('dimensionless')
    
    # Estimate gamma (Cp/Cv) assuming Cv ~ Cp/1.25 if no direct Cv is available
    specific_heat_v = specific_heat_p / 1.25
    gamma = specific_heat_p / specific_heat_v
    
    # Coefficient of thermal expansion (1/°R) as an approximation
    alpha = Q_(1 / (temperature.magnitude * 60), 1 / unitReg.degR)
    
    # Thermal diffusivity (ft²/s)
    thermal_diffusivity = thermal_conductivity / (density * specific_heat_p) * unitReg.foot**2 / unitReg.second

    # Surface Tension (lbf/ft)
    SurfaceTens = Q_(fluid.SurfAtTdegR(temperature.magnitude), unitReg.pound_force / unitReg.foot)
    
    # Return all properties in English Engineering units
    return (viscosity, specific_heat_p, gamma, thermal_conductivity, density, prandtl, alpha, thermal_diffusivity, SurfaceTens)

def get_film_properties(name, temperature_R):
    # Ensure inputs are in the correct units for RocketProps methods
    temperature = Q_(temperature_R.magnitude, unitReg.degR)
    
    # Create an fluid object
    fluid = get_prop(name)

    # Latent heat of vaporization
    h_fg = Q_(fluid.HvapAtTdegR(temperature.magnitude), unitReg.BTU / unitReg.pound)
    
    # Saturation Properties
    Psat = Q_(fluid.PvapAtTdegR(temperature.magnitude), unitReg.psi)
    Tsat = Q_(fluid.TdegRAtPsat(Psat.magnitude), unitReg.degR)
    return (h_fg, Psat, Tsat)

def DarcyFrictionFactor(reynoldsNum, surfaceRoughness, hydroDiameter):
    if reynoldsNum < 2000: # laminar
        return 64/reynoldsNum
    
    # if reynoldsNum < 4000:
        # print("Warning: Flow is transitional, Darcy friction factor may not be accurate.")
    # turbulent
    f_func = lambda f: -2*np.log10(surfaceRoughness/hydroDiameter/3.7 + 2.51/(reynoldsNum*np.sqrt(f))) - 1/np.sqrt(f)
    f = fsolve(f_func, 0.05)    # Darcy friction factor
    try:
        return f[0]
    except IndexError:
        return f
    
def FrictionPressureLoss(f: float, length: float, hydroDiameter: float, density: float, velocity: float) -> float:
    return length/hydroDiameter * f * density * velocity**2 / 2