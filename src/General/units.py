import pint
from ambiance import Atmosphere

# -------------- Unit Shit -------------- #
unitReg = pint.UnitRegistry()
unitReg.default_system = 'US'
unitReg.formatter.default_format = "~P"  # Compact unit formatting
unitReg.define("lbmol = 453.59237 * mol")  # 1 lbmol = 453.59237 mol (since 1 lb = 453.59237 g)

Q_ = unitReg.Quantity

# constants
R_UNIVERSAL = Q_(10.731577089016, unitReg.psi * unitReg.foot**3 / (unitReg.lbmol * unitReg.degR))  # Universal gas constant in ft·lbf/(lbmol·°R)
PRESCOTT_ALT = Q_(5400, unitReg.feet)
prescottAtm = Atmosphere(PRESCOTT_ALT.to(unitReg.meter).magnitude)
PRESCOTT_PRESSURE = Q_(prescottAtm.pressure[0], unitReg.pascal)
PRESCOTT_TEMP = Q_(70, unitReg.degF)