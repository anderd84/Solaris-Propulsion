import Nozzle.rao as rao
import numpy as np
import matplotlib.pyplot as plt
import matrix_viewer as mv
from fluids.gas import PrandtlMeyerFunction, SpHeatRatio, MachAngle
import fluids.gas as gas
import Nozzle.plots as nozplt
import Nozzle.nozzle as nozzle
import Nozzle.analysis as analysis

gamma = SpHeatRatio(1.17)
Me = 2.942
Te = np.deg2rad(-3.5)
Mt = 1
PbPc = 0.015


Pc = 300
PambPc = 6.75/Pc

print(np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2]))


data = rao.GenerateInputMatrix(np.linspace(2.5,3.5,10), np.deg2rad(np.linspace(-6, -1, 10)), SpHeatRatio(1.17), .002)
rao.GenerateInputChart(data)

plt.show()

# mv.view(data)
# mv.show()








Md = rao.CalculateMachD(Me, Te, gamma, PbPc)
areaRatio, Cf, lengthRatio = rao.CalculatePlugMetrics(Me, Te, Md, gamma, steps=100000)

# areaRatio = 30.7
# lengthRatio = 1.9554

print(areaRatio, Cf, lengthRatio)
# input("Press Enter to continue...")


Tt = rao.CalculateThroatAngle(Me, Te, Mt, gamma)
phit = Tt  + (np.pi/2 - MachAngle(Mt))
# length = 1.9554
# expansionRatio = 5

length = lengthRatio
expansionRatio = areaRatio

controlSurface: np.ndarray[rao.CharacteristicPoint] = rao.GetControlSurfaceProperties(Me, Te, length, gamma, 50)
expansionFan: np.ndarray[rao.CharacteristicPoint] = rao.GenerateExpansionFan(Me, Mt, Tt, gamma , 50)

field = rao.GenerateFlowField(expansionFan, controlSurface, gamma)

Rt = np.sqrt(1 - (1/expansionRatio*np.cos(Tt)))

print(np.pi/np.cos(Tt)*(1-Rt**2))

print(f"Throat Angle: {Tt * 180/np.pi - 90}, Throat Radius: {Rt}")

cont = rao.CalculateContour(field, Rt, Tt)

# controlSurface: np.ndarray[rao.CharacteristicPoint] = rao.GetControlSurfaceProperties(Me, Te, length, gamma, 100)
# expansionFan: np.ndarray[rao.CharacteristicPoint] = rao.GenerateExpansionFan(Me, Mt, Tt, gamma , 200)

# field = rao.GenerateFlowField(expansionFan, controlSurface, gamma)

field = rao.PruneField(field)

# field = rao.PruneUnderContour(field, cont)

x = np.array([[p.x for p in row] for row in field])
r = np.array([[p.r for p in row] for row in field])
mach = np.array([[p.mach for p in row] for row in field])
theta = np.array([[p.theta for p in row] for row in field])
alpha = np.array([[p.alpha for p in row] for row in field])

# mv.view(np.nan_to_num(x), "X")

# mv.view(np.nan_to_num(r), "R")

# mv.view(np.nan_to_num(alpha), "alpha")

# np.nan_to_num(mach, copy=False)
# mv.view(mach, "Mach")

# np.nan_to_num(theta, copy=False)
# mv.view(np.rad2deg(theta), "Theta")

# mv.view(np.transpose(np.array([cx, cy])), "Contour")

# mv.show()

# nozplt.show3d(cont)
nozplt.WriteContourTXT(cont, "contour.txt")
fieldplot = nozplt.CreateNonDimPlot()
fieldplot = nozplt.PlotField(fieldplot, field)
fieldplot = nozplt.PlotContour(fieldplot, cont, Rt, Tt)
# fieldplot.axes[0].legend(['Control Surface', 'Throat', 'Nozzle'], prop={'size': 12})
# fieldplot.axes[0].set_xlim(cont[-1].x, cont[0].x)
# nozplt.LiveContour(field, Rt, Tt, fieldplot)
# internal,external = nozzle.InternalPreExpansion(Me, Te, Mt, Tt, gamma, (cont[-1].x, cont[-1].r), expansionRatio)
# fieldplot.axes[0].plot(internal[0,:], internal[1, :], 'r')
# fieldplot.axes[0].plot(external[0,:], external[1, :], 'b')
# fieldplot.axes[0].plot(internal, external, '.')

print(f"Cf: {Cf}")
FT = Cf*Pc


# Rt = np.sqrt(1 - (1/expansionRatio*np.cos(phit)))
# fieldplot.axes[0].plot([0, (1 - Rt)*np.tan(phit)], [1, Rt], '-g', linewidth=2) # Throat
plt.show()



contour2 = nozzle.RaoContourFormat(cont)
# analysis.CalculateSimpleField(contour2, 0.005, 0, gamma, Mt, Tt, steps=100)

mat, stream = analysis.CalculateComplexField(contour2, PambPc, PbPc, gas.Gas(gamma+0, 287), Mt, Tt, 25, 0, 2)
fieldplot2 = nozplt.CreateNonDimPlot()
fieldplot2 = nozplt.PlotContour(fieldplot2, cont, Rt, Tt)
# ic(stream)
fieldplot2.axes[0].plot([p.x for p in stream], [p.r for p in stream], '--b', linewidth=1.5)
analysis.PlotCharacteristicLines(fieldplot2, mat)
Me = np.sqrt((PambPc**(-1/gamma[5]) - 1)/gamma[2])
thetaExit = Tt + gas.PrandtlMeyerFunction(Me, gamma) - gas.PrandtlMeyerFunction(Mt, gamma)
# fieldGrid = analysis.GridifyComplexField(mat, np.array([]))
fieldplot2.axes[0].legend(['Throat', 'Nozzle', 'Free Boundary'], loc=1, prop={'size': 12})
fieldplot2.axes[0].grid('on', linestyle='--')
# analysis.PlotFieldData(fieldplot2, fieldGrid)
# fieldplot2.axes[0].plot([0, 2*np.cos(thetaExit)], [1, 1+2*np.sin(thetaExit)], 'g')
plt.show()