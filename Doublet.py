from scipy.optimize import fsolve
from dataclasses import dataclass
import matplotlib.pyplot as plt
from pint import UnitRegistry
from scipy import integrate
from scipy import linalg
import numpy as np
import scipy as sp 
rg = np.random.default_rng()
ureg = UnitRegistry()
ureg.default_system = 'US'
Q_ = ureg.Quantity
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
        Velocity = Cd * np.sqrt(2 * Pressure_Diff/ self.rho )
        return Velocity.to(ureg.feet / ureg.second)
    
    
    def Area(self, Cd, Pressure_Diff) -> float:
        """Calculates area from Dr. Whites Eq 16.2
        Args:
            Cd (float): Coefficient of Discharge
            Pressure_Diff (float): Pressure difference through injector (Estimated metric) in psi
        Returns:
            float: Total Orifice Area needed for Propellant
        """        
        Area = (self.mdot / (Cd * np.sqrt(2*self.rho* Pressure_Diff )))
        return Area.to(ureg.inch**2)
    
    
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
        Number = np.round(Tot_Area/Hole_Area, decimals=0)
        return Number
    
#Constants
CD_drill = 0.7 #Constant for Sharp Edged Orifices
g0 = Q_(32.174, ureg.foot / ureg.second**2)
#Design Input
#gamma defined as angle after faceplate off the axial direction. Positive defined as away from center of Aerospike
mdots = np.array([5.29, 2.21]) #LOX_CORE, FUEL_CORE
Film_Cooling = np.array([0.08, 0.08]) #Outer Film Cooling, Inner Film Cooling
Exit_Angle = 15 #degrees. Design Constraint
Pressure_Drop_Fuel = 0.2
Pressure_Drop_Lox = 0.2
Pressure_Chamber = Q_(300, ureg.force_pound / ureg.inch**2)
Doublet_Diameter_LOX = Q_(0.125, ureg.inch)
Doublet_Diameter_Fuel = Q_(0.0625, ureg.inch) #Guess 
#Equations and Setting classes
OX_CORE = PROP(gamma=0, mdot=mdots[0], rho=56.794)
FUEL_CORE = PROP(gamma = 0, mdot = mdots[1], rho=51.15666)
def func(gamma_FUEL):
    return -1 * np.tan(np.deg2rad(Exit_Angle)) + \
    (OX_CORE.mdot.magnitude * OX_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Lox).magnitude * np.sin(np.deg2rad(OX_CORE.gamma.magnitude)) \
    + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel).magnitude * np.sin(gamma_FUEL))/ \
    (OX_CORE.mdot.magnitude * OX_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Lox).magnitude * np.cos(np.deg2rad(OX_CORE.gamma.magnitude)) \
    + FUEL_CORE.mdot.magnitude * FUEL_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel).magnitude * np.cos(gamma_FUEL))    
print(OX_CORE.mdot * OX_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Lox)  )
print(FUEL_CORE.mdot * FUEL_CORE.Velocity(CD_drill, Pressure_Chamber * Pressure_Drop_Fuel))  
gamma_FUEL = fsolve(func, np.deg2rad(5)) 
FUEL_CORE.gamma = Q_(np.rad2deg(gamma_FUEL[0]), ureg.degrees) 
print(f"Fuel Angle is {FUEL_CORE.gamma:.3f}")
OUT_FILM_C = PROP(gamma = 30, mdot = Film_Cooling[0]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
IN_FILM_C = PROP(gamma = -30, mdot = Film_Cooling[1]* FUEL_CORE.mdot, rho = FUEL_CORE.rho)
OX_CORE_Holes = OX_CORE.Number(Doublet_Diameter_LOX,CD_drill, Pressure_Chamber * (Pressure_Drop_Lox))
print(OX_CORE_Holes)
FUEL_CORE_Diameter = np.sqrt((FUEL_CORE.Area(CD_drill, Pressure_Chamber * (Pressure_Drop_Fuel)) / OX_CORE_Holes) * 4 / np.pi)
print(FUEL_CORE_Diameter) #Currently gives us roughly 0.082 inches. Closest Drill size is 2.08 mm diamter which is a Drill size of 45
FUEL_CORE_Holes = FUEL_CORE.Number(Q_(2.08, ureg.mm),CD_drill, Pressure_Chamber * (Pressure_Drop_Fuel))
print(FUEL_CORE_Holes)
print(f"Total FUEL Doublets Velocity: {FUEL_CORE.Velocity(CD_drill, Pressure_Chamber * (Pressure_Drop_Fuel)):.3f~}")
print(f"Total FUEL Orifice Area Doublets: {FUEL_CORE.Area(CD_drill, Pressure_Chamber * (Pressure_Drop_Fuel)):.3f~}")
print(f"Total OX Doublets Velocity: {OX_CORE.Velocity(CD_drill, Pressure_Chamber * (Pressure_Drop_Lox)):.3f~}")
print(f"Total OX Orifice Area Doublets: {OX_CORE.Area(CD_drill, Pressure_Chamber * (Pressure_Drop_Lox)):.3f~}")
FUEL_mdot_tot = FUEL_CORE.mdot + OUT_FILM_C.mdot + IN_FILM_C.mdot
print(f"For Total Fuel mdot in normal units: {FUEL_mdot_tot:.3f~}")
print(f"For Total Fuel mdot in commie units: {FUEL_mdot_tot.to(ureg.kilogram / ureg.second):.3f~}")


# PLOTTING SHIT BELOW


# Constants and parameters
Points = 1000
di = 6.5
ri = di / 2
Spacing = 0.5  # in
Rgamma_lox = 1.5  # in
gamma_lox = 15  # degrees
gamma_fuel = -12  # degrees
Chamber_Cowl_r = 0.5  # in
Past_Peak = 1.15

#These following constants are only used for the shitty Aerospike contour I made for proof of concept for my impingement formula
r1 = 3.50
r2 = 2.50
angle1 = 41.81
angle2 = 83.62

# Initial positions
startX, startY = 0, 1
startX1, startY1 = 2, 1  # Adjusted based on given code
BaseX = np.linspace(startX, startX1, Points)
BaseY = np.linspace(startY, startY1, Points)

# Calculating end points and start points for the arcs
endX1 = startX1 + r1 * np.sin(np.radians(angle1))
endY1 = startY1 + r1 * (1 - np.cos(np.radians(angle1)))
startX2, startY2 = endX1, endY1
endX2 = startX2 + r2 * (np.cos(np.radians(90 - angle2 / 2)) - np.cos(np.radians(90 + angle2 / 2)))
endY2 = startY2 + r2 * (np.sin(np.radians(90 - angle2 / 2)) - np.sin(np.radians(90 + angle2 / 2)))
startX3, startY3 = endX2, endY2
endX3 = startX3 + r1 * (np.cos(3*np.pi/2) - np.cos(np.radians(270 - angle1)))
endY3 = startY3 + r1 * (np.sin(3*np.pi/2) - np.sin(np.radians(270 - angle1)))

# Defining theta ranges for the arcs
thetaRange1 = np.linspace(0, np.radians(angle1), Points)
thetaRange2 = np.linspace(np.radians(90 + angle2 / 2), np.radians(90 - angle2 / 2), Points)
thetaRange3 = np.linspace(np.radians(angle1), 0, Points)

# Defining arcs
arc1_x = startX1 + r1 * np.sin(thetaRange1)
arc1_y = startY1 + r1 * (1 - np.cos(thetaRange1))
arc2_x = startX2 + r2 * (np.cos(np.radians(90 - angle2 / 2)) + np.cos(thetaRange2))
arc2_y = startY2 + r2 * (-np.sin(np.radians(90 - angle2 / 2)) + np.sin(thetaRange2))
arc3_x = endX3 + r1 * np.cos(3 * np.pi / 2 - thetaRange3)
arc3_y = endY3 + r1 * (1 + np.sin(3 * np.pi / 2 - thetaRange3))

# Combining all segments
x_profile = np.concatenate([BaseX, arc1_x, arc2_x, arc3_x])
y_profile = np.concatenate([BaseY, arc1_y, arc2_y, arc3_y])

# Plotting the aerospike nozzle contour
plt.figure()
plt.plot(x_profile, y_profile, linewidth=2)
plt.axhline(0, linewidth=2, color="b")
plt.axhline(Rgamma_lox, linestyle="--")
plt.axhline(Rgamma_lox + Spacing, linestyle="--")
Peaky, Peak_Point = max(y_profile), np.argmax(y_profile)
Peakx = x_profile[Peak_Point]
x_graph = np.linspace(0, max(x_profile), Points)
gamma_lox_line = np.tan(np.radians(gamma_lox)) * x_graph + Rgamma_lox
plt.plot(x_graph, gamma_lox_line, "--")
gamma_fuel_line = np.tan(np.radians(gamma_fuel)) * x_graph + Rgamma_lox + Spacing
plt.plot(x_graph, gamma_fuel_line, "--")

# Calculations for impingement point and chamber contour
gamma_fuel = abs(gamma_fuel)
x = Spacing * np.sin(np.radians(90 - gamma_fuel)) / np.sin(np.radians(gamma_lox + gamma_fuel)) * np.sin(np.radians(90 - gamma_lox))
y = Spacing * np.sin(np.radians(90 - gamma_fuel)) / np.sin(np.radians(gamma_lox + gamma_fuel)) * np.cos(np.radians(90 - gamma_lox)) + Rgamma_lox
plt.plot(x, y, "o")

ChamberX = x_graph / x_graph[-1] * Peakx * Past_Peak
ChamberY = np.ones(len(ChamberX)) * ri
thetaRange4 = np.linspace(np.radians(90), 0, Points)
ChamberArcX = ChamberX[-1] + Chamber_Cowl_r * np.cos(thetaRange4)
ChamberArcY = ChamberY[-1] + Chamber_Cowl_r * (-1 + np.sin(thetaRange4))
Chamber_ContourX = np.concatenate([ChamberX, ChamberArcX])
Chamber_ContourY = np.concatenate([ChamberY, ChamberArcY])
plt.plot(Chamber_ContourX, Chamber_ContourY, "k", linewidth=2)

# Additional plotting for resultant line
yprime = (ri + Peaky) / 2
xprime = Peakx
plt.plot(xprime, yprime, "o")
tan_resultant = (yprime - y) / (xprime - x)
resultant_y_intercept = y - tan_resultant * x
ResultantX = x_graph / x_graph[-1] * (xprime - x) * Past_Peak + x
ResultantY = tan_resultant * ResultantX + resultant_y_intercept
plt.plot(ResultantX, ResultantY, "g", linewidth=2)

# Finishing touches
plt.legend(['Spike Contour', 'Centerline', 'Gamma_{OX} Straight Line', 'Gamma_{FUEL} Straight Line', 'Gamma_{OX} Angled Line', 'Gamma_{FUEL} Angled Line',
            'Impingement Point', 'Chamber Contour', 'Aim Point', 'Resultant Line'], loc="best")
plt.xlabel('Distance Along Engine Axis (inches)')
plt.ylabel('Radius (inches)')
plt.axis('equal')
plt.title('Side View Contour of an Aerospike Nozzle')
plt.grid(True)
plt.show()
