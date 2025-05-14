from alive_progress import alive_bar
from cooling import material
from cooling.material import CoolantType, DomainMaterial, MaterialType
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from cooling import domain
from nozzle import plots
import general.design as DESIGN
from nozzle import plug
from nozzle import analysis
from general.units import Q_, unitReg

TEST_INLET = Q_(522, unitReg.psi)  # Test inlet pressure
TEST_INJECTOR = Q_(170, unitReg.psi)  # Test injector pressure

Re = Q_(3.2, unitReg.inch)
exhaust = DESIGN.exhaustGas
print(exhaust.stagTemp)
print(exhaust.stagPress)
print(DESIGN.Fuel_Total)
print(DESIGN.Oxidizer_Total)

cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
Rt = outputData["radiusThroat"]
Tt = outputData["thetaThroat"]
Re = outputData["radiusLip"]

overchoke = plug.getOverchokeDist(Re, Rt, Tt, DESIGN.chokePercent)

plugC, straightLength, plugCoolL, plugCoolU = plug.GenerateDimPlug(cont, Rt, Tt, Re, Q_(6.3, unitReg.inch), Q_(1.5, unitReg.inch))
cowlC, cowlCoolL, cowlCoolU = plug.GenerateDimCowl(Rt, Tt, Re, straightLength, DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke)
chamberC, aimpoint = plug.GenerateDimChamber(Rt, Tt, Re, Q_(6.3, unitReg.inch), DESIGN.chamberInternalRadius, DESIGN.wallThickness, overchoke, Q_(1.5, unitReg.inch))

highmesh = domain.DomainMC(-7.5, 3.3, 8.5, 1.5, .01)

highmesh.DefineMaterials(cowlC, chamberC, plugC, 10)

outerloop = highmesh.NewCoolantLoop(Q_(.025, 'inch'), 300, DESIGN.Fuel_Total, CoolantType.H2O)
highmesh.AssignCoolantFlow(domain.CoolingChannel(cowlCoolU, cowlCoolL), False, Q_(TEST_INJECTOR, unitReg.psi), outerloop)

surfaceRough = fsolve(lambda eps: (highmesh.ChannelPressureDrop(Q_(eps, 'inch'), outerloop) - (TEST_INLET - TEST_INJECTOR)).m, .1)
print(f"Surface roughness: {surfaceRough[0]} inch")