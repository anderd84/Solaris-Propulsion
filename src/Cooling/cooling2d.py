from matplotlib import pyplot as plt, patches
import numpy as np
import os
from scipy.optimize import fsolve
from icecream import ic

#from Doublet import OX_CORE, FUEL_CORE
import General.design as DESIGN
from fluids.fluid import CD_drill, Pressure_Drop_Fuel, \
                         Pressure_Drop_Lox,  pm, get_fluid_properties, CP
from Injector.Doublet_Functions import spike_contour
from Injector.InjectorCad import injector_cad_write
from Injector.Drill import drill_approximation
from General.units import Q_, unitReg

#First step always is to update doublet.py file and run beforehand to grab all mdot and density values at injector side

# -------------- Design Inputs -------------- #
NumberofChannels = 30
ChannelShape =  np.array([.25, .25]) #inches Width,Height


def wall_thermal_conductivity(T):
    # Returns thermal conductivity of GR-COP-42 based on temperature
    T = T.to(unitReg.degK)
    return -9E-05*T^2 + 0.091*T + 327.88    # W/m/K

# Code for generating convection coefficients

def combustion_convection(mu, c_p, Pr, P_0, c_star, D_star, A_star, A, r_c,
                          Ma, T_wg, T_0g, gamma):
    # Bartz's Correlation for combustion side
    # Inputs (Combustion gases)
    mu = mu.to(unitReg.pound / unitReg.foot / unitReg.second)  # Dynamic viscosity at stagnation conditions
    c_p = c_p.to(unitReg.BTU / unitReg.pound / unitReg.degR)    # Specific heat at stagnation conditions
    # Pr = 2.45       # Prandtl number at stagnation conditions
    P_0 = P_0.to(unitReg.pound / unitReg.inch**2)   # Stagnation pressure
    c_star = c_star.to(unitReg.feet / unitReg.s) # Characteristic velocity
    D_star = D_star.to(unitReg.feet)      # Throat diameter (as hydraulic diameter, 4*A_c/P)
    A_star = A_star.to(unitReg.inch**2)      # Throat area
    A = A.to(unitReg.inch**2)           # Nozzle area at location of interest
    r_c = r_c.to(unitReg.inch)   # Throat radius of curvature
    # Ma = 0          # Local Mach number
    T_wg = T_wg.to(unitReg.degR)     # Hot side wall temperature
    T_0g = T_0g.to(unitReg.degR)  # Hot gas stagnation temperature
    omega = 0.6     # for diatomic gases
    # gamma = 1.145   # Ratio of specific heats, assumed to be constant

    # Calculations
    sigma = 1 / ((1 / 2 * T_wg / T_0g * (1 + (gamma - 1) / 2 * Ma**2) + 1 / 2)**(0.8 - 0.2 * omega) * (1 + (gamma - 1) / 2 * Ma**2)**(0.2 * omega))
    return 0.026 / D_star**0.2 * mu**0.2 / Pr**0.6 * c_p * (P_0 / c_star)**0.8 * (D_star / r_c)**0.1 * (A_star / A)**0.9 * sigma  # Convective heat transfer coefficient

def f_equation(f, *data):
    epsilon, D_h, Re_D = data
    return -2*np.log10(epsilon/D_h/3.7 + 2.51/(Re_D*np.sqrt(f))) - 1/np.sqrt(f)

# TO DO: design file epsilon, pull properties using temperature as a function argument
def internal_flow_convection(epsilon, m_dot_c, NumberofChannels, rho, v, A_c, P, Pr, mu, k_c, mu_s):
    # Gnielinski/Sieder & Tate for channel side
    # Inputs (Cooling channel)
    epsilon = epsilon.to(unitReg.inch)     # Surface roughness
    m_dot_c = m_dot_c.to(unitReg.pound / unitReg.second) / NumberofChannels     # Coolant mass flow rate through one channel
    rho = rho.to(unitReg.pound / unitReg.inch**3)  # Coolant density
    v = v.to(unitReg.feet / unitReg.s)      # Coolant velocity along channel
    D_h = D_h.to(unitReg.inch)     # Hydraulic diameter of cooling channel
    # Pr = 1      # Prandtl number (likely an array)
    mu = mu.to(unitReg.pound / unitReg.foot / unitReg.second)   # Dynamic viscosity (likely an array)
    k_c = k_c.to(unitReg.BTU / unitReg.foot / unitReg.hour / unitReg.degR)    # Thermal conductivity of coolant
    mu_s = mu_s.to(unitReg.pound / unitReg.foot / unitReg.second)    # Dynamic viscosity at the heat transfer boundary surface temperature

    # Calculations
    D_h = 4*A_c/P
    Re_D = 4*m_dot_c/np.pi/D_h/mu   # Reynold's number
    # Re_D = Re_D.to(unitReg.inch / unitReg.inch)
    # Darcy Friction Factor
    if Re_D < 2300:
        f = 64/Re_D # Laminar flow
    else:
        data = (epsilon, D_h, Re_D) # Arguments for fsolve
        f = fsolve(f_equation, 0.05, args=data)  # Turbulent flow
    # Convection coefficient (Gnielinski)
    if Pr >= 0.7 and Pr <= 16700 and Re_D >= 3000 and Re_D <= 5000000:   # Check that properties fit restrictions for Gnielinski
        Nu_D = f/8*(Re_D - 1000)*Pr / (1 + 12.7*(f/8)**0.5 * (Pr**(2/3) - 1))   # Nusselt number of fully developed flow
    else:    # Use Sieder & Tate otherwise
        Nu_D = 0.027*Re_D**0.8*Pr**(1/3)*(mu/mu_s)**0.14    # Nusselt number of fully developed flow
    return Nu_D*k_c/D_h      # Convective heat transfer coefficient

def free_convection(beta, T_s, T_infinity, P_atm, D_outer):
    # Properties
    beta = beta.to(1 / unitReg.degR)
    T_s = T_s.to(unitReg.degR)
    T_infinity = T_infinity.to(unitReg.degR)
    P_atm = P_atm.to(unitReg.pound / unitReg.inch**2)
    D_outer = D_outer.to(unitReg.inch)
    g = Q_(32.174, unitReg.feet / unitReg.s**2)
    (viscosity, _, _, k_f, density, Pr, _, _, _) = get_fluid_properties('Air', T_infinity, P_atm)
    nu = viscosity/density  # Dynamic viscosity, air property lookup
    # Calculations
    Gr_D = Q_(g*beta*(T_s - T_infinity)*D_outer**3/nu^2, 
              unitReg.inch / unitReg.inch)
    Gr_D.to_base_units()
    Ra_D = Q_(Gr_D*Pr, '')
    if Ra_D < 10^12:
        Nu_D = (0.60 + 0.387*Ra_D^(1/6)/(1 + (0.559/Pr)^(9/16))^(8/27))^2
    else:
        Nu_D = -999999999999
        print("Raleigh number exceeds restriction")
    return Nu_D*k_f/D_outer

def lambda_equation(Lambda, *data):
    Re_g = data
    return 1.930*np.log10(Re_g*np.sqrt(Lambda)) - 1/np.sqrt(Lambda)

def film_cooling(m_dot_g, m_dot_c, u_g, u_c, P_cc, D_cc, c_p_g, mu_g, Pr_g, 
                 rho_g, M_g, sigma_g, mu_c, c_c_l, h_fg, T_c_1, T_c_sat, rho_c_l, 
                 M_c, h_g, D_c):
    # Inputs
    # viscosity, specific_heat_p, gamma, thermal_conductivity, density, prandtl, alpha, thermal_diffusivity, SurfaceTens
    mu_c, c_c_l, gamma, _, rho_c_l, _, _, _, _ = get_fluid_properties(DESIGN.fuelName, temperature_R, pressure_psi)
    mu_g, c_p_g, gamma, _, density, Pr_g, alpha, _, sigma_g = get_fluid_properties(DESIGN.oxName, temperature_R, pressure_psi)
    # m_dot_g   Combustion gas mass flow rate
    # m_dot_c   Film coolant mass flow rate
    # u_g       Combustion gas velocity
    # u_c       Film coolant velocity
    # P_cc      Chamber pressure
    # D_cc      Chamber diameter (might be hydraulic?)
    # c_p_g     Combustion gas specific heat at constant pressure
    # mu_g      Combustion gas dynamic viscosity
    # Pr_g      Combustion gas Prandtl number
    # rho_g     Combustion gas density
    # M_g       Combustion gas molecular weight
    # sigma_g   Combustion gas surface tension
    # mu_c      Film coolant dynamic viscosity
    # c_c_l     Film coolant specific heat as liquid
    # h_fg      Film coolant enthalpy of vaporization
    # T_c_1     Film coolant initial temperature
    # T_c_sat   Film coolant saturation temperature
    # rho_c_l   Film coolant liquid density
    # M_c       Film coolant molecular weight
    # h_g       Combustion gas convection coefficient
    # D_c       Film cooling channel diameter
    G_g = rho_g*u_g # Combustion gas mass velocity
    G_mean = G_g*(u_g - u_c)/u_g    # Mean mass velocity
    Re_g = G_mean*D_cc/mu_g # Combustion gas Reynolds number
    data = Re_g
    Lambda = fsolve(lambda_equation, 0.1, args=data) # Friction factor
    e_t = 0.1   # Must be found from testing data?
    K_t = 1 + 4*e_t # Corrective turbulence factor
    St = Lambda/2/(1.20 + 11.8*np.sqrt(Lambda/2)*(Pr_g - 1)*(Pr_g)^(-1/3))  # Stanton number
    h_o = G_mean*c_p_g*St*K_t # Dry wall convection coefficient
    h_fg_star = h_fg + (T_c_sat - T_c_1)*c_c_l
    epsilon = 0 # Must take from Leckner's data for spectral radiative properties of combustion gases
    sigma = 5.67*10^-8  # Stefan-Boltzmann constant
    Q_dot_rad = sigma*epsilon*np.pi/4*D_cc**2*(T_g**4 - T_c**4) # TO DO: Might need to change area
    q_dot_rad = Q_dot_rad/(np.pi/4*D_cc**2)
    q_dot_conv = h_g*(T_g-T_c)
    Q_dot_conv = q_dot_conv # TO DO
    m_dot_v = (q_dot_conv + q_dot_rad)/h_fg_star    # Coolant evaporation rate per area
    F = m_dot_g/(rho_g*u_g) # Blowing ratio
    St_o = St/(np.log(1 + F/St*(M_g/M_c)^0.6)/(F/St*(M_g/M_c)^0.6)) # Transpiration-corrected Stanton number
    h = St_o*rho_c_l*u_c*c_c_l  # Convection coefficient
    Re_c = rho_c_l*u_c*D_c/mu_c # Film coolant Reynolds number
    a = 2.31*10^-4*Re_c^-0.35
    Re_cfilm = 250*np.log(Re_c) - 1265
    E_m = 1- Re_cfilm/Re_c  # Maximum entrainment fraction
    We = rho_g*u_g^2*D_cc/sigma_g*((rho_c_l-rho_g)/rho_g)^(0.25)   # Weber number TO DO: Check for correct variables (coolant vs combustion gases)
    E = E_m*np.tanh(a*We^1.25)  # Entrainment fraction
    Gamma_c = m_dot_c*(1-E)/(np.pi*D_cc)    # TO DO: Check that this is right
    L_c = Gamma_c/m_dot_v   # Liquid film cooled length
    # Q_dot_conv = h
# coolMesh.mesh.z
# Testing
# """ temperature_R = 700 * unitReg.degR  # Temperature in Rankine (~80°F)
# pressure_psi = 14.7 * unitReg.psi  # Pressure in psi (1 atmosphere)
# properties = get_fluid_properties('n-Dodecane', temperature_R.magnitude, 
#                                   pressure_psi.magnitude)
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