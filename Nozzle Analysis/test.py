import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import rao
import gas
import matrix_viewer as mv

Me = 3
Te = np.deg2rad(-8)
gamma = gas.SpHeatRatio(1.23)
PbPc = 0

# print(rao.InputDataGenerate2(Me, Te, gamma, PbPc))
print(rao.CalculatePlugMetrics(Me, Te, rao.CalculateMachD(Me, Te, gamma, PbPc), gamma, steps=100))

data = rao.GenerateInputMatrix(np.linspace(2.1, 4.4, 25), np.deg2rad(np.linspace(-8, 0, 9)), gamma, PbPc)

rao.GenerateInputChart(data)
