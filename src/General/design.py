from rocketcea.cea_obj import CEA_Obj
from ambiance import Atmosphere

from General.units import Q_, unitReg
from General.units import R_UNIVERSAL
from fluids.gas import Gas

# eninge inputs

# props
oxName = "LOX"
fuelName = "RP-1"
OFratio = 2

# geometry
chamberPressure = Q_(300, unitReg.psi)
totalmdot = Q_(8, unitReg.pound / unitReg.sec)
maxRadius = Q_(3, unitReg.inch)
wallThickness = Q_(0.25, unitReg.inch)

# injector inputs
chamberInternalRadius = maxRadius - wallThickness #Internal Diameter of Chamber
Spacing = Q_(0.55, unitReg.inch)  #Spacing between centear of impingement Holes
oxHoleRadius = Q_(1.25, unitReg.inch)  #Radial distance between centerline and LOX hole
percentFilmCooling = 0.15 #Outer Film Cooling Percentage
filmCoolingSpacing = Q_(0.60, unitReg.inch) #inches Outer
oxDoubletDiameter = Q_(0.0625, unitReg.inch)  #Design choise for DOublet Diameter size (Need to look more properly into it as 1/4 holes might make vaporization time too long)\
oxImpingeAngle = Q_(25, unitReg.degrees)
filmImpingeAngle = Q_(25, unitReg.degrees)
oxDewarPressure = Q_(22, unitReg.psi)
prescottAmbientTemp = Q_(70, unitReg.degF)
prescottAmbientPressure = Q_(12.2, unitReg.force_pound / unitReg.inch**2)

# Nozzle inputs
designAltitude = Q_(20000, unitReg.feet)
designAtm = Atmosphere(designAltitude.to(unitReg.meter).magnitude)

designAmbientPressure = Q_(designAtm.pressure[0], unitReg.pascal)
lengthMax = Q_(5.5, unitReg.inch)
basePressure = Q_(20, unitReg.psi)

# chamber derived
Combustion=CEA_Obj(oxName=oxName, fuelName=fuelName);
chamberTemp = Q_(Combustion.get_Tcomb(Pc=chamberPressure.magnitude, MR=OFratio), unitReg.degR)

# throat derived
_, _, _, molarWeightThroat, gamma = Combustion.get_IvacCstrTc_ThtMwGam(Pc=chamberPressure.magnitude, MR=OFratio, eps=1)
molarWeightThroat = Q_(molarWeightThroat, unitReg.pound / unitReg.lbmol)
R_throat = (R_UNIVERSAL / molarWeightThroat).to(unitReg.foot * unitReg.pound_force / (unitReg.pound * unitReg.degR))

# exhaust derived
exhaustGas = Gas(gamma, R_throat, P0=chamberPressure, T0=chamberTemp)
chokeArea = exhaustGas.getChokedArea(totalmdot).to(unitReg.inch**2)

# nozzle design table
plugDesignTable = {"throatArcRadFactor": .1, "convergeAngle": 30, "turnArcRadFactor": 2, "straightAngle": 10, "lipAngle": 15}