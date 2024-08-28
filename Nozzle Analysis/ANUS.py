import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from pint import UnitRegistry

import gas
import angelino as Angelino
import rao as Rao

ureg = UnitRegistry()
ureg.default_system = 'US'

# DESIGN INPUTS
gamma = 1.4
pressureChamber = 300 * ureg.psi
machExit = 4.5
radiusCowlLip = 3 * ureg.inch
areaThroat = .3 * ureg.inch**2
radiusBase = 0 * ureg.inch


areaRatio, Cf, lengthRatio, aplot, rplot, lplot, tplot, mplot =  Rao.InputDataGenerate(2.4, np.deg2rad(-8.25), 1.23, .001, 0)
# areaRatio, Cf, lengthRatio, aplot, rplot, lplot, tplot, mplot = Rao.InputDataGenerate(2.4, -1, 1.4, .00001, .01)

print()
print(areaRatio, Cf, lengthRatio)
print(mplot[-1], np.rad2deg(tplot[-1]))

plt.plot(mplot, tplot)
plt.grid(True)
plt.show()

# plt.plot(rplot, aplot) # v
# plt.grid(True)
# plt.show()

# plt.plot(rplot, np.rad2deg(tplot)) # ^

# plt.ylim((-22, -6))
# plt.xlim((0, 1))
# plt.gca().invert_xaxis()
# plt.gca().set_xticks(np.linspace(0, 1, 11))
# plt.grid(True, 'both')
# plt.show()


# Rao.GenerateInputChart(np.linspace(2,4.4,25), np.linspace(-8, 0, 9), 1.4, .001, .0001)