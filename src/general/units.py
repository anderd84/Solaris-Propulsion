import pint
from ambiance import Atmosphere
from enum import IntEnum

# -------------- Unit Shit -------------- #
unitReg = pint.application_registry.get()
unitReg.default_system = 'US'
unitReg.formatter.default_format = "~P"  # Compact unit formatting
unitReg.define("lbmol = 453.59237 * mol")  # 1 lbmol = 453.59237 mol (since 1 lb = 453.59237 g)

Q_ = unitReg.Quantity

# constants
R_UNIVERSAL = Q_(10.731577089016, unitReg.psi * unitReg.foot**3 / (unitReg.lbmol * unitReg.degR))  # Universal gas constant in ft·lbf/(lbmol·°R)
G0 = Q_(32.174, unitReg.foot/unitReg.second**2)
PRESCOTT_ALT = Q_(5400, unitReg.feet)
prescottAtm = Atmosphere(PRESCOTT_ALT.to(unitReg.meter).magnitude)
PRESCOTT_PRESSURE = Q_(prescottAtm.pressure[0], unitReg.pascal)
PRESCOTT_TEMP = Q_(70, unitReg.degF)

# ez units
R = unitReg.degR
FT = unitReg.feet
IN = unitReg.inch
PSI = unitReg.psi
LBF = unitReg.pound_force
LBM = unitReg.pound
S = unitReg.seconds

# directions
class Direction(IntEnum):
    LEFT = 0
    UPPER = 1
    UP = 1
    LOWER = 2
    DOWN = 2
    RIGHT = 3