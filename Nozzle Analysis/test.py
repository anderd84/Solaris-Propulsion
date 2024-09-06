import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import rao
import gas

Me = 2.4
Te = np.deg2rad(-8.25)
gamma = gas.SpHeatRatio(1.23)
PbPc = 0

# print(rao.InputDataGenerate2(Me, Te, gamma, PbPc))
print(rao.CalculatePlugMetrics(Me, Te, gamma, PbPc, steps=100000))
plt.show()