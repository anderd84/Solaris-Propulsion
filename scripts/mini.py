import numpy as np
from scipy.optimize import fsolve
from icecream import ic

from fluids.fluid import get_fluid_properties, get_film_properties
from general.units import Q_, R_UNIVERSAL, unitReg
import general.design as DESIGN
from general.design import exhaustGas

epsilon = DESIGN.epsilon
fuelname = DESIGN.fuelName
pressure_stagnation = DESIGN.chamberPressure
temperature_stagnation = DESIGN.chamberTemp
specific_heat_stagnation = DESIGN.heat_capacity_chamber
viscosity_stagnation = DESIGN.viscosity_chamber 
Prandtl_stagnation = DESIGN.Prandtl_chamber
Combustion = DESIGN.Combustion
cstar = DESIGN.c_star
R_gas = DESIGN.R_throat #TODO gotta make this better so it's not constant R assumption
gamma_throat = DESIGN.gamma #TODO gotta make this better so it's not constant R assumption
A_star = DESIGN.chokeArea
mdot_tot = DESIGN.totalmdot
Fuel_Total = DESIGN.Fuel_Total

# Area calculation for internal external nozzle to flash to gas at throat
T_throat = Q_(680, unitReg.degR)
P = Q_(120, unitReg.psi)
(_, Psat, _) = get_film_properties(DESIGN.fuelName, T_throat)
(_, _, gamma, _, rho, _, _, _, _) = get_fluid_properties(DESIGN.fuelName, P, T_throat)
q = Q_(2, unitReg.pound/unitReg.second)/rho
A = Q_(np.pi/4*(0.8**2 - 0.76**2), unitReg.inch**2)
P_1 = P
P_2 = Psat # Q_(5.5, unitReg.psi)
A_star = (A*q*np.sqrt(rho/(2*A**2*P_1 - 2*A**2*P_2 + q**2*rho))).to(unitReg.inch**2)
r_star = np.sqrt(A_star/np.pi).to(unitReg.inch)
 
m_dot = Q_(2, unitReg.pound/unitReg.second)
A_star = np.sqrt((1/A**2 - 2*(P_2-P_1)*rho/m_dot**2)**-1).to(unitReg.inch**2)
r_star = np.sqrt(A_star/np.pi).to(unitReg.inch)
 
# Area calculation for internal external nozzle for Mach 1 at throat
(_, _, gamma, _, rho, _, _, _, _) = get_fluid_properties(DESIGN.fuelName, Psat, T_throat)
v_throat = (m_dot/(rho*A_star)).to(unitReg.feet/unitReg.second)
P_e = Q_(5.5, unitReg.psi)
M = Q_(175, unitReg.gram/unitReg.mole).to(unitReg.pound/unitReg.lbmol)
R_g = R_UNIVERSAL/M
rho = Psat/R_g/T_throat
P_0 = Psat + 0.5*rho*v_throat**2
a = np.sqrt(gamma*R_g*T_throat).to(unitReg.feet/unitReg.second)
Ma = v_throat/a
A = A_star
A_star = A*((gamma + 1)/2)**((gamma + 1)/(2*(gamma - 1)))*Ma/(1 + (gamma - 1)/2*Ma**2)**((gamma + 1)/(2*(gamma - 1)))
r_star = np.sqrt(A_star/np.pi).to(unitReg.inch)
Ma_e = np.sqrt(((P_0/P_e)**((gamma - 1)/gamma) - 1) * 2/(gamma - 1))
A_e = A_star*((gamma + 1)/2)**(-(gamma + 1)/(2*(gamma - 1)))*(1 + (gamma - 1)/2*Ma_e**2)**((gamma + 1)/(2*(gamma - 1)))/Ma_e

ic(A_star)
ic(A_e)