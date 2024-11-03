import matplotlib.pyplot as plt

from Cooling import domain
from Nozzle import plots

coolmesh = domain.DomainMC.LoadFile("coolmesh.pkl")

fig = plots.CreateNonDimPlot()
coolmesh.ShowMaterialPlot(fig)

plt.show()