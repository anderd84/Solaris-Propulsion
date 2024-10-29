from matplotlib import pyplot as plt, patches
import numpy as np
import os
from scipy.optimize import fsolve
from icecream import ic

#from Doublet import OX_CORE, FUEL_CORE
from fluids.fluid import Q_, unitreg, CD_drill, Pressure_Drop_Fuel, \
                         Pressure_Drop_Lox,  pm, get_fluid_properties, CP
from Injector.Doublet_Functions import spike_contour
from Injector.InjectorCad import injector_cad_write
from Injector.Drill import drill_approximation
from generalInput import Q_, unitreg

#First step always is to update doublet.py file and run beforehand to grab all mdot and density values at injector side

# -------------- Design Inputs -------------- #
NumberofChannels = 30
ChannelShape =  np.array([.25, .25]) #inches Width,Height


# Code for generating convection coefficients

def combustion_convection(mu, c_p, Pr, P_0, c_star, D_star, A_star, A, r_c,
                          Ma, T_wg, T_0g, gamma):
    # Bartz's Correlation for combustion side
    # Inputs (Combustion gases)
    mu = mu.to(unitreg.pound / unitreg.foot / unitreg.second)  # Dynamic viscosity at stagnation conditions
    c_p = c_p.to(unitreg.BTU / unitreg.pound / unitreg.degR)    # Specific heat at stagnation conditions
    # Pr = 2.45       # Prandtl number at stagnation conditions
    P_0 = P_0.to(unitreg.pound / unitreg.inch**2)   # Stagnation pressure
    c_star = c_star.to(unitreg.feet / unitreg.s) # Characteristic velocity
    D_star = D_star.to(unitreg.feet)      # Throat diameter (as hydraulic diameter, 4*A_c/P)
    A_star = A_star.to(unitreg.inch**2)      # Throat area
    A = A.to(unitreg.inch**2)           # Nozzle area at location of interest
    r_c = r_c.to(unitreg.inch)   # Throat radius of curvature
    # Ma = 0          # Local Mach number
    T_wg = T_wg.to(unitreg.degR)     # Hot side wall temperature
    T_0g = T_0g.to(unitreg.degR)  # Hot gas stagnation temperature
    omega = 0.6     # for diatomic gases
    # gamma = 1.145   # Ratio of specific heats, assumed to be constant

    # Calculations
    sigma = 1 / ((1 / 2 * T_wg / T_0g * (1 + (gamma - 1) / 2 * Ma**2) + 1 / 2)**(0.8 - 0.2 * omega) * (1 + (gamma - 1) / 2 * Ma**2)**(0.2 * omega))
    return 0.026 / D_star**0.2 * mu**0.2 / Pr**0.6 * c_p * (P_0 / c_star)**0.8 * (D_star / r_c)**0.1 * (A_star / A)**0.9 * sigma  # Convective heat transfer coefficient

def func(f, *data):
    epsilon, D_h, Re_D = data
    return -2*np.log10(epsilon/D_h/3.7 + 2.51/(Re_D*np.sqrt(f))) - 1/np.sqrt(f)

def internal_flow_convection(epsilon, m_dot_c, NumberofChannels, rho, v, A_c, P, Pr, mu, k_c, mu_s):
    # Gnielinski/Sieder & Tate for channel side
    # Inputs (Cooling channel)
    epsilon = epsilon.to(unitreg.inch)     # Surface roughness
    m_dot_c = m_dot_c.to(unitreg.pound / unitreg.second) / NumberofChannels     # Coolant mass flow rate through one channel
    rho = rho.to(unitreg.pound / unitreg.inch**3)  # Coolant density
    v = v.to(unitreg.feet / unitreg.s)      # Coolant velocity along channel
    D_h = D_h.to(unitreg.inch)     # Hydraulic diameter of cooling channel
    # Pr = 1      # Prandtl number (likely an array)
    mu = mu.to(unitreg.pound / unitreg.foot / unitreg.second)   # Dynamic viscosity (likely an array)
    k_c = k_c.to(unitreg.BTU / unitreg.foot / unitreg.hour / unitreg.degR)    # Thermal conductivity of coolant
    mu_s = mu_s.to(unitreg.pound / unitreg.foot / unitreg.second)    # Dynamic viscosity at the heat transfer boundary surface temperature

    # Calculations
    D_h = 4*A_c/P
    Re_D = 4*m_dot_c/np.pi/D_h/mu   # Reynold's number
    # Re_D = Re_D.to(unitreg.inch / unitreg.inch)
    # Darcy Friction Factor
    if Re_D < 2300:
        f = 64/Re_D # Laminar flow
    else:
        data = (epsilon, D_h, Re_D) # Arguments for fsolve
        f = fsolve(func, 0.05, args=data)  # Turbulent flow
    # Convection coefficient (Gnielinski)
    if Pr >= 0.7 and Pr <= 16700 and Re_D >= 3000 and Re_D <= 5000000:   # Check that properties fit restrictions for Gnielinski
        Nu_D = f/8*(Re_D - 1000)*Pr / (1 + 12.7*(f/8)**0.5 * (Pr**(2/3) - 1))   # Nusselt number of fully developed flow
    else:    # Use Sieder & Tate otherwise
        Nu_D = 0.027*Re_D**0.8*Pr**(1/3)*(mu/mu_s)**0.14    # Nusselt number of fully developed flow
    return Nu_D*k_c/D_h      # Convective heat transfer coefficient

def free_convection(beta, T_s, T_infinity, P_atm, D_outer):
    # Properties
    beta = beta.to(1 / unitreg.degR)
    T_s = T_s.to(unitreg.degR)
    T_infinity = T_infinity.to(unitreg.degR)
    P_atm = P_atm.to(unitreg.pound / unitreg.inch**2)
    D_outer = D_outer.to(unitreg.inch)
    g = Q_(32.174, unitreg.feet / unitreg.s**2)
    (viscosity, _, _, k_f, density, Pr, _, _, _) = get_fluid_properties('Air', T_infinity, P_atm)
    nu = viscosity/density  # Dynamic viscosity, air property lookup
    # Calculations
    Gr_D = Q_(g*beta*(T_s - T_infinity)*D_outer**3/nu^2, unitreg.inch / unitreg.inch)
    Gr_D.to_base_units()
    Ra_D = Q_(Gr_D*Pr, '')
    if Ra_D < 10^12:
        Nu_D = (0.60 + 0.387*Ra_D^(1/6)/(1 + (0.559/Pr)^(9/16))^(8/27))^2
    else:
        Nu_D = -999999999999
        print("Raleigh number exceeds restriction")
    return Nu_D*k_f/D_outer

# Testing
# """ temperature_R = 700 * unitreg.degR  # Temperature in Rankine (~80°F)
# pressure_psi = 14.7 * unitreg.psi  # Pressure in psi (1 atmosphere)
# properties = get_fluid_properties('n-Dodecane', temperature_R.magnitude, pressure_psi.magnitude)
# h_g = combustion_convection()
# #fluids = CP.get_global_param_string('fluids_list')

# # Output the results
# (viscosity, specific_heat_p, gamma, thermal_conductivity, density, prandtl, 
#  quality, phase, alpha_SI, thermal_diffusivity) = properties

# print(f"Viscosity: {viscosity:.6f}")
# print(f"Specific Heat (Cp): {specific_heat_p:.6f}")
# print(f"Gamma (Cp/Cv): {gamma:.6f}")
# print(f"Thermal Conductivity: {thermal_conductivity:.6f}")
# print(f"Density: {density:.6f}")
# print(f"Prandtl Number: {prandtl:.6f}")
# print(f"Quality: {quality:.6f}")
# print(f"Phase: {phase:.6f}")
# print(f"Coefficient of Thermal Expansion (1/K): {alpha_SI:.6f}")
# print(f"Thermal Diffusivity: {thermal_diffusivity:.6f}") """