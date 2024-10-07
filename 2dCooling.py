from matplotlib import pyplot as plt, patches
from Doublet_Functions import spike_contour
from InjectorCad import injector_cad_write
from Drill import drill_approximation
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from fluid import Q_, ureg, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox,  pm, get_fluid_properties, CP
#from Doublet import OX_CORE, FUEL_CORE
import numpy as np
import os

#First step always is to update doublet.py file and run before hand to grab all mdot and density values at injector side

# -------------- Design Inputs -------------- #
NumberofChannels = 30
ChannelShape =  np.array([.25, .25]) #inches Width,Height




# Code for generating convection coefficients

# Bartz's Correlation for combustion side
# Inputs (Combustion gases)
mu = 0.0000905  # Dynamic viscosity at stagnation conditions
c_p = 5488.9    # Specific heat at stagnation conditions
Pr = 2.45       # Prandtl number at stagnation conditions
P_0 = 2054600   # Stagnation pressure
c_star = 1789.1 # Characteristic velocity
D_star = 2      # Throat diameter (as hydraulic diameter, 4*A_c/P)
A_star = 2      # Throat area
A = 2           # Nozzle area at location of interest
r_c = 0.00254   # Throat radius of curvature
Ma = 0          # Local Mach number
T_wg = 2900     # Hot side wall temperature
T_0g = 3328.15  # Hot gas stagnation temperature
omega = 0.6     # for diatomic gases
gamma = 1.145   # Ratio of specific heats, assumed to be constant

# Calculations
sigma = 1 / ((1 / 2 * T_wg / T_0g * (1 + (gamma - 1) / 2 * Ma**2) + 1 / 2)**(0.8 - 0.2 * omega) *
              (1 + (gamma - 1) / 2 * Ma**2)**(0.2 * omega))
h_g = 0.026 / D_star**0.2 * mu**0.2 / Pr**0.6 * c_p * (P_0 / c_star)**0.8 * (D_star / r_c)**0.1 * (A_star / A)**0.9 * sigma  # Convective heat transfer coefficient

print(f"Combustion side convective heat transfer coefficient: {h_g:.2f} W/m^2/K")


# Gniliesnski/Sieder & Tate for channel side
# Inputs (Cooling channel)
f = 0.000   # Friction factor
rho = Q_(1000, ureg.pound / ureg.foot**3)  # Coolant density
v = Q_(10, ureg.ft / ureg.second**2)      # Coolant velocity along channel
A_h = Q_(10, ureg.ft**2) 
P_h = Q_(10, ureg.ft ) 
D_h = 4*A_h/P_h     # Hydraulic diameter of cooling channel
Pr = 1      # Prandtl number (likely an array)
mu = 0.01   # Dynamic viscosity (likely an array)
k_c = 20    # Thermal conductivity of coolant
mu_s = 1    # Dynamic viscosity at the heat tranfer boundary surface temperature

# Calculations
Re_D = rho*v*D_h/mu     # Reynold's number
# Re_D = 4*m_dot_c/numpy.pi/D_h/mu      # Alternate formula
if Pr >= 0.7 and Pr <= 16700 and Re_D >= 3000 and Re_D <= 5000000:   # Check that properties fit restrictions for Gnielinski
    Nu_D = f/8*(Re_D - 1000)*Pr / (1 + 12.7*(f/8)**0.5 * (Pr**(2/3) - 1))   # Nusselt number of fully developed flow
    h_c = Nu_D*k_c/D_h      # Convective heat transfer coefficient
else:    # Use Sieder & Tate otherwise
    Nu_D = 0.027*Re_D**0.8*Pr**(1/3)*(mu/mu_s)**0.14 # Nusselt number of fully developed flow
    h_c = Nu_D*k_c/D_h      # Convective heat transfer coefficient

# TEsting SHit
#fluids = CP.get_global_param_string('fluids_list')
#print(f"Available fluids in CoolProp: {fluids}")



temperature_R = 700 * ureg.degR  # Temperature in Rankine (~80°F)
pressure_psi = 14.7 * ureg.psi  # Pressure in psi (1 atmosphere)

# Call the function to get properties
properties = get_fluid_properties('n-Dodecane',temperature_R.magnitude, pressure_psi.magnitude)

# Output the results
(viscosity, specific_heat_p, gamma, thermal_conductivity, density, prandtl, 
 quality, phase, alpha_SI, thermal_diffusivity) = properties

print(f"Viscosity: {viscosity:.6f}")
print(f"Specific Heat (Cp): {specific_heat_p:.6f}")
print(f"Gamma (Cp/Cv): {gamma:.6f}")
print(f"Thermal Conductivity: {thermal_conductivity:.6f}")
print(f"Density: {density:.6f}")
print(f"Prandtl Number: {prandtl:.6f}")
print(f"Quality: {quality:.6f}")
print(f"Phase: {phase:.6f}")
print(f"Coefficient of Thermal Expansion (1/K): {alpha_SI:.6f}")
print(f"Thermal Diffusivity: {thermal_diffusivity:.6f}")