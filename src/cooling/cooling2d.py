import numpy as np
from scipy.optimize import fsolve

from fluids.fluid import get_fluid_properties
from general.units import Q_, unitReg
import general.design as DESIGN
from general.design import exhaustGas

epsilon = DESIGN.epsilon
fuelname = DESIGN.fuelName
Combustion = DESIGN.Combustion
pressure_stagnation = DESIGN.chamberPressure
temperature_stagnation = DESIGN.chamberTemp
specific_heat_stagnation = DESIGN.heat_capacity_chamber
viscosity_stagnation = DESIGN.viscosity_chamber 
Prandtl_stagnation = DESIGN.Prandtl_chamber
cstar = DESIGN.c_star
R_gas = DESIGN.R_throat #TODO gotta make this better so it's not constant R assumption
gamma_throat = DESIGN.gamma #TODO gotta make this better so it's not constant R assumption
A_star = DESIGN.chokeArea
NumberofChannels = DESIGN.NumberofChannels
mdot_tot = DESIGN.totalmdot
Fuel_Total = DESIGN.Fuel_Total

#First step always is to update doublet.py file and run beforehand to grab all mdot and density values at injector side

# -------------- Design Inputs -------------- #

ChannelShape =  np.array([.25, .25]) #inches Width,Height

mdotperchannel = Fuel_Total / NumberofChannels





# @unitReg.wraps((unitReg.degR, unitReg.psi), (unitReg.degR, unitReg.degR, None, unitReg.psi, unitReg.psi, unitReg.inch**2, unitReg.inch, unitReg.inch))
def heatcoolant(Tprev, Tcell, resSet, Pcell):
    TRsum, Rsum = resSet.getSums()
    (_, cp, _, _, _, _, _, _, _) = get_fluid_properties(fuelname, Tcell, Pcell)

    Tprev = Tprev.to(unitReg.degR)

    Tnew = ((TRsum + Tprev*mdotperchannel*cp)/(mdotperchannel*cp + Rsum)).to(unitReg.degR)

    return Tnew

def depresscoolant(Tcell, Pprev, Pcell, channelArea, channelHydroD, DeltaL):
    (mu, _, _, _, rho, _, _, _, _) = get_fluid_properties(fuelname, Tcell, Pcell)
    A_c = channelArea # Cooling channel cross-sectional area, height should be variable in the future
    D_h = channelHydroD   # Hydraulic diameter
    Re_D = mdotperchannel*D_h/mu/A_c    # Reynolds number
    DeltaL = Q_(DeltaL, unitReg.inch)   # Step size

    f_func = lambda f: -2*np.log10(epsilon/D_h/3.7 + 2.51/(Re_D*np.sqrt(f))) - 1/np.sqrt(f)
    f = fsolve(f_func, 0.05)    # Darcy friction factor
    DeltaP = (f*DeltaL/D_h*rho*(mdotperchannel/rho/A_c)**2/2).to(unitReg.psi)   # Pressure drop

    Pcell = Pcell.to(unitReg.psi)
    Pprev = Pprev.to(unitReg.psi)

    Pnew = Pprev - DeltaP

    return Pnew












def conduction_rp1(Temp):
    Temp = Q_(Temp.magnitude, unitReg.degR)
    (_, _, _, thermal_conductivity, _,_, _, _, _) =get_fluid_properties(fuelname, Temp.to(unitReg.degR), pressure_stagnation.to(unitReg.psi))
    return thermal_conductivity



def conduction_grcop(Temp):
    Temp = Temp.to(unitReg.degK)
    if Temp.magnitude>1000:
        Temp = Q_(1000, unitReg.degK)
    conduction_coefficient = Q_(-9E-05*Temp.magnitude**2 + 0.091*Temp.magnitude + 327.88 , unitReg.watt / (unitReg.meter * unitReg.degK))    # W/m/K
    return conduction_coefficient.to(unitReg.BTU / unitReg.foot / unitReg.hour / unitReg.degR)


# Code for generating convection coefficients

def calculate_nozzle_area(A_star, M, gamma):
    # Calculate area A using the rearranged equation
    A = A_star * ((gamma + 1) / 2) ** (- (gamma + 1) / (2 * (gamma - 1))) * (1 + (gamma - 1) / 2 * M**2) ** ((gamma + 1) / (2 * (gamma - 1))) * (1 / M)
    return A



def combustion_convection(Node_Temp, Velocity):
    Node_Temp = Q_(Node_Temp.magnitude, unitReg.degR)
    Velocity = Q_(Velocity.magnitude, unitReg.foot / unitReg.second)
    gamma = exhaustGas.SimpleHarmonicGamma(Node_Temp).g
    Temp = temperature_stagnation.to(unitReg.degR)
    sonic_velocity = np.sqrt((gamma * R_gas * Node_Temp)).to(unitReg.foot / unitReg.second)
    Mach = Velocity / sonic_velocity

    # Define constants and convert variables with units
    mu = viscosity_stagnation.to(unitReg.pound / (unitReg.foot * unitReg.second))  # Dynamic viscosity at stagnation conditions
    c_p = specific_heat_stagnation.to(unitReg.BTU / (unitReg.pound * unitReg.degR))  # Specific heat at stagnation conditions
    Pr = Prandtl_stagnation  # Prandtl number (dimensionless)
    P_0 = pressure_stagnation.to(unitReg.pound_force / unitReg.ft**2)  # Convert pressure to lbf/ft²
    c_star = cstar.to(unitReg.foot / unitReg.second)  # Characteristic velocity in ft/s
    A_star = DESIGN.chokeArea.to(unitReg.ft**2)  # Convert throat area to ft²

    #TODO make this section use the actual R_e and R_t
    P_T = 2 * np.pi * (design.R_T + design.R_E)  # Perimeter in feet
    D_star = (4 * A_star / P_T).to(unitReg.ft)  # Throat diameter as hydraulic diameter in feet

    # Calculate the nozzle area at the location of interest and convert to ft²
    A = calculate_nozzle_area(A_star, Mach, gamma).to(unitReg.ft**2)
    r_c = Q_(0.1, unitReg.inch).to(unitReg.ft)  # Throat radius of curvature in feet#TODO make it way more accurate

    # Local wall and gas stagnation temperatures
    T_wg = Node_Temp.to(unitReg.degR)  # Hot side wall temperature
    T_0g = Temp  # Hot gas stagnation temperature
    omega = 0.6  # for diatomic gases

    # Calculate sigma correction factor
    sigma = 1 / ((0.5 * T_wg / T_0g * (1 + (gamma - 1) / 2 * Mach**2) + 0.5)**(0.8 - 0.2 * omega) *
                 (1 + (gamma - 1) / 2 * Mach**2)**(0.2 * omega))

    g_c = Q_(32.174, unitReg.pound * unitReg.ft / (unitReg.pound_force * unitReg.s**2))  # Gravity constant to convert lbf to lbm

    h_g = (0.026 / D_star**0.2 * mu**0.2 / Pr**0.6 * c_p * (P_0 / (c_star * g_c))**0.8 *
           (D_star / r_c)**0.1 * (A_star / A)**0.9 * sigma)

    # Return h_g in the desired heat transfer coefficient units
    return h_g.to(unitReg.BTU / (unitReg.foot**2) / unitReg.hour / unitReg.degR)




#Velocity_test = Q_(3000, unitReg.foot / unitReg.second)
#Temp_test = Q_(5800, unitReg.degR)
#ic(combustion_convection2(Temp_test,Velocity_test))

def internal_flow_convection(Node_Temp, Node_Pressure, channelArea, channelHydroD):
    # Gnielinski/Sieder & Tate for channel side
    # Inputs (Cooling channel)
    Temp = Q_(Node_Temp.magnitude, unitReg.degR)    # Not sure about this being the right temperature for properties
    Pressure = Q_(Node_Pressure.magnitude, unitReg.psi)
    (mu, _,_, k_c, _, Pr, _, _, _) = get_fluid_properties(fuelname, Temp, Pressure) # Coolant property lookup
    
    
    
     
    # epsilon = epsilon.to(unitReg.inch)     # Surface roughness from Design Table
    m_dot_c = mdotperchannel.to(unitReg.pound / unitReg.second)     # Coolant mass flow rate through one channel from Design Table
    mu = mu.to(unitReg.pound / unitReg.foot / unitReg.second)   # Dynamic viscosity #! From RP-1 Grab values from ROcket Prop
    k_c = k_c.to(unitReg.BTU / unitReg.foot / unitReg.hour / unitReg.degR)    # Thermal conductivity of coolant #! From RP-1 Grab values from ROcket Prop
    mu_s = mu.to(unitReg.pound / unitReg.foot / unitReg.second)    # Dynamic viscosity at the heat transfer boundary surface temperature NEEDS CORRECTING

    # Calculations
    A_c = channelArea # Cooling channel cross-sectional area, height should be variable in the future
    D_h = channelHydroD   # Hydraulic diameter
    Re_D = (m_dot_c/A_c)*D_h/mu   # Reynolds number
    # Re_D = Re_D.to(unitReg.inch / unitReg.inch)
    # Darcy Friction Factor
    if Re_D < 2300:
        f = 64/Re_D # Laminar flow
    else:
        f_func = lambda f: -2*np.log10(epsilon/D_h/3.7 + 2.51/(Re_D*np.sqrt(f))) - 1/np.sqrt(f)
        f = fsolve(f_func, 0.05)  # Turbulent flow
    # Convection coefficient (Gnielinski)
    f = f[0]
    if Pr >= 0.5 and Pr <= 2000 and Re_D >= 3000 and Re_D <= 5000000:   # Check that properties fit restrictions for Gnielinski
        Nu_D = f/8*(Re_D - 1000)*Pr / (1 + 12.7*(f/8)**0.5 * (Pr**(2/3) - 1))   # Nusselt number of fully developed flow
    elif Pr >= 0.7 and Pr <= 16700 and Re_D >= 10000:    # Use Sieder & Tate otherwise
        print("Using Sieder-Tate, not ready yet")
        Nu_D = 0.027*Re_D**0.8*Pr**(1/3)*(mu/mu_s)**0.14    # Nusselt number of fully developed flow
    else:
        print("Pr/Re is wack")
        Nu_D = 0.027*Re_D**0.8*Pr**(1/3)*(mu/mu_s)**0.14    # Use Sieder-Tate anyway
    return (Nu_D*k_c/D_h).to((unitReg.BTU / unitReg.foot**2 / unitReg.hour / unitReg.degR))      # Convective heat transfer coefficient

def plug_convection_coefficient(P, v, x):
    P = P.to(unitReg.psi)   # Local exhaust gas pressure
    v = v.to(unitReg.feet / unitReg.second) # Local exhaust gas velocity
    x = x.to(unitReg.feet)  # Arc length along plug from throat
    (_, mu, k, Pr) = Combustion.get_Exit_Transport(P, DESIGN.OFratio, pressure_stagnation/P, 0, 0)  # Exit viscosity, thermal conductivity, Prandtl number
    (_, _, rho) = Combustion.get_Densities(P, DESIGN.OFratio, pressure_stagnation/P, 0, 0)  # Exit density
    mu = mu.to(unitReg.pound / unitReg.foot / unitReg.second)   # Viscosity
    k = k.to(unitReg.BTU / unitReg.foot / unitReg.hour / unitReg.degR)  # Thermal conductivity
    x = x.to(unitReg.feet)  # Arc length along plug
    rho = rho.to(unitReg.pound / unitReg.foot**3)   # Combustion gas density
    L = (4*DESIGN.chokeArea/(2*np.pi*(DESIGN.R_E + DESIGN.R_T)) + x).to(unitReg.feet)   # Throat hydraulic diameter + arc length
    Re_x = rho*v*L/mu   # Reynolds number
    C_f = 0.455/(np.log(0.06*Re_x))**2  # Friction coefficient
    h = k/L*C_f/2*Re_x*Pr/(1 + 12.7*(Pr**(2/3) - 1)*np.sqrt(C_f/2)) # Convection coefficient
    # T_ref = 
    # rho_ref = 
    # mu_ref = 
    # h = 0.023*(rho_ref/rho)**0.8*(mu_ref/mu)**0.2*Re_x**(-0.2)/Pr**0.6/(Re_x*Pr)  # Stanton-Reynolds number relation for turbulent flow
    return h.to(unitReg.BTU / (unitReg.foot**2) / unitReg.hour / unitReg.degR)

def free_convection(T_s, T_infinity, P_atm, D_outer):
    # Properties
    beta = 1/T_s.to(1 / unitReg.degR)
    T_s = T_s.to(unitReg.degR)
    T_infinity = T_infinity.to(unitReg.degR)
    P_atm = P_atm.to(unitReg.pound / unitReg.inch**2)
    D_outer = D_outer.to(unitReg.inch)
    g = Q_(32.174, unitReg.foot / unitReg.s**2)
    (viscosity, _, _, k_f, density, Pr, _, _, _) = get_fluid_properties('Air', T_infinity, P_atm)
    nu = viscosity/density  # Dynamic viscosity, air property lookup
    # Calculations
    Gr_D = Q_(g*beta*(T_s - T_infinity)*D_outer**3/nu**2, 
              unitReg.inch / unitReg.inch)
    Gr_D.to_base_units()
    Ra_D = Q_(Gr_D*Pr, '')
    if Ra_D < 10**12:
        Nu_D = (0.60 + 0.387*Ra_D**(1/6)/(1 + (0.559/Pr)**(9/16))**(8/27))**2
    else:
        Nu_D = -999999999999
        print("Raleigh number exceeds restriction")
    return Nu_D*k_f/D_outer

def get_film_length(m_dot_c, A_c, P, u_g, u_c, rho_g, mu_c, c_c_l, h_fg,
                    T_g, T_c_1, T_c_sat, rho_c_l, sigma_g, h_g, D_c):
    # Inputs
    m_dot_c = m_dot_c.to(unitReg.pound / unitReg.second)    # Combustion gas mass flow rate
    A_c = A_c.to (unitReg.feet**2)  # Chamber cross-sectional area
    P = P.to (unitReg.feet) # Chamber wetted perimeter
    u_g = u_g.to(unitReg.feet / unitReg.second) # Combustion gas velocity
    u_c = u_c.to(unitReg.feet / unitReg.second) # Film coolant velocity
    rho_g = rho_g.to(unitReg.pound/ unitReg.foot**3)    # Combustion gas density
    mu_c = mu_c.to(unitReg.pound / unitReg.foot / unitReg.second)   # Film coolant dynamic viscosity
    c_c_l = c_c_l.to(unitReg.BTU / unitReg.pound / unitReg.degR)    # Film coolant specific heat as liquid
    h_fg = h_fg.to(unitReg.BTU / unitReg.pound / unitReg.degR)  # Film coolant enthalpy of vaporization
    T_g = T_g.to(unitReg.degR)  # Combustion gas temperature
    T_c_1 = T_c_1.to(unitReg.degR)  # Film coolant initial temperature
    T_c_sat = T_c_sat.to(unitReg.degR)  # Film coolant saturation temperature
    rho_c_l = rho_c_l.to(unitReg.pound/ unitReg.foot**3)    # Film coolant liquid density
    sigma_g = sigma_g.to(unitReg.pound / unitReg.foot)  # Combustion gas surface tension
    h_g = h_g.to((unitReg.BTU / unitReg.foot**2 / unitReg.hour / unitReg.degR)) # Combustion gas convection coefficient
    D_c = D_c.to(unitReg.foot)  # Film cooling channel diameter

    D_h = 4*A_c/P
    h_fg_star = h_fg + (T_c_sat - T_c_1)*c_c_l  # Heat of vaporization + enthalpy change from subcooled to saturated
    epsilon = 0 # Combustion gas emissivity Must take from Leckner's data for spectral radiative properties of combustion gases
    sigma = Q_(1.713441*10**-9, unitReg.BTU / unitReg.hour / unitReg.foot**2 / unitReg.degR**4)  # Stefan-Boltzmann constant
    q_dot_rad = sigma*epsilon*(T_g**4 - T_c_1**4) # Radiative heat flux, assuming negligible right now
    q_dot_conv = h_g*(T_g - T_c_1)  # Convective heat flux
    q_dot_tot = q_dot_rad + q_dot_conv  # Total heat rate into coolant
    m_dot_v = q_dot_tot/h_fg_star   # Evaporated coolant per unit area
    Re_c = rho_c_l*u_c*D_c/mu_c # Film coolant Reynolds number
    a = 2.31*10**-4*Re_c**-0.35
    Re_cfilm = 250*np.log(Re_c) - 1265
    E_m = 1- Re_cfilm/Re_c
    We = rho_g*u_g**2*D_h/sigma_g   # Weber number TO DO: Check for correct variables, especially D and sigma
    E = E_m*np.tanh(a*We**1.25)
    Gamma_c = m_dot_c*(1-E)/(np.pi*2*DESIGN.chamberInternalRadius)    # TO DO: Check that this is right
    L_c = Gamma_c/m_dot_v
    return L_c

def film_cooling(m_dot_g, A_c, P, u_g, u_c, mu_g, Pr_g, 
                 rho_g, M_g, c_c_l, rho_c_l, M_c):
    # Inputs
    m_dot_g = m_dot_g.to(unitReg.pound / unitReg.second)    # Combustion gas mass flow rate
    A_c = A_c.to (unitReg.feet**2)  # Chamber cross-sectional area
    P = P.to (unitReg.feet) # Chamber wetted perimeter
    u_g = u_g.to(unitReg.feet / unitReg.second) # Combustion gas velocity
    u_c = u_c.to(unitReg.feet / unitReg.second) # Film coolant velocity
    mu_g = mu_g.to(unitReg.pound / unitReg.foot / unitReg.second)   # Combustion gas dynamic viscosity
    # Pr_g      Combustion gas Prandtl number
    rho_g = rho_g.to(unitReg.pound/ unitReg.foot**3)    # Combustion gas density
    M_g = M_g.to(unitReg.pound / unitReg.lbmol)    # Combustion gas molecular weight
    c_c_l = c_c_l.to(unitReg.BTU / unitReg.pound / unitReg.degR)    # Film coolant specific heat as liquid
    rho_c_l = rho_c_l.to(unitReg.pound/ unitReg.foot**3)    # Film coolant liquid density
    M_c = M_c.to(unitReg.pound / unitReg.lbmol)    # Film coolant molecular weight

    D_h = 4*A_c/P
    G_g = rho_g*u_g # Combustion gas mass velocity
    G_mean = G_g*(u_g - u_c)/u_g    # Mean mass velocity
    Re_g = G_mean*D_h/mu_g # Combustion gas Reynolds number
    lambda_func = lambda Lambda: 1.930*np.log10(Re_g*np.sqrt(Lambda)) - 0.537 - 1/np.sqrt(Lambda)
    Lambda = fsolve(lambda_func, 0.1) # Friction factor
    # e_t = 0.1   # RMS turbulence fraction Must be found from testing data?
    # K_t = 1 + 4*e_t # Corrective turbulence factor
    St = Lambda/2/(1.20 + 11.8*np.sqrt(Lambda/2)*(Pr_g - 1)*(Pr_g)**(-1/3))  # Stanton number
    # h_o = G_mean*c_p_g*St*K_t # Dry wall convection coefficient, not used since we choose Bartz
    F = m_dot_g/(rho_g*u_g) # Blowing ratio
    St_o = St/(np.log(1 + F/St*(M_g/M_c)**0.6)/(F/St*(M_g/M_c)**0.6)) # Transpiration-corrected Stanton number
    h = St_o*rho_c_l*u_c*c_c_l  # Convection coefficient
    return h
