import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from pint import UnitRegistry

import gas
import nozzle

ureg = UnitRegistry()
ureg.default_system = 'US'

# DESIGN INPUTS
gamma = 1.4
pressureChamber = 300 * ureg.psi
machExit = 4.5
radiusCowlLip = 3 * ureg.inch
areaThroat = .3 * ureg.inch**2
radiusBase = 0 * ureg.inch


areaRatio, Cf, lengthRatio, machA, thetaA, mplot, tplot = nozzle.RaoInputDataGenerate(2.4, -8.25, 1.23, .001, 0)

print(areaRatio, Cf, lengthRatio, machA, thetaA)

