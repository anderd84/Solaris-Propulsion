import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import rao
import gas
import matrix_viewer as mv

mach = np.linspace(1, 500, 100)
gamma = gas.SpHeatRatio(1.4)

nus = gas.PrandtlMeyerFunction(mach, gamma)
plt.plot(mach, np.rad2deg(nus))
plt.grid(True)
plt.show()