from re import T
import rao
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matrix_viewer as mv

gamma = 1.23
Me = 2.4
Te = np.deg2rad(-8.25)
Mt = 1.3
Tt = rao.CalculateThroatAngle(Me, Te, Mt, gamma)
PbPc = 0
length = 1.164
expansionRatio = 3.81



controlSurface: np.ndarray[rao.CharacteristicPoint] = rao.GetControlSurfaceProperties(Me, Te, length, gamma, 10)
expansionFan: np.ndarray[rao.CharacteristicPoint] = rao.GenerateExpansionFan(Me, Te, Mt, Tt, gamma , 10)

csx = [p.x for p in controlSurface]
csr = [p.r for p in controlSurface]
cst = [p.theta for p in controlSurface]
csm = [p.mach for p in controlSurface]

# print()
# print(cst, csm)
# # plt.plot(cst)
# plt.plot(csx, csm)

efx = [p.x for p in expansionFan]
efr = [p.r for p in expansionFan]
eft = [p.theta for p in expansionFan]
efm = [p.mach for p in expansionFan]

# print()
# print(eft, efm)
# # plt.plot(eft)
# plt.plot(np.rad2deg(eft))

# plt.scatter(csx, csr)
# plt.scatter(efx, efr)
# plt.grid(True)
# plt.show()

field = rao.GenerateFlowField(expansionFan, controlSurface, gamma)


x = np.array([[p.x for p in row] for row in field])
r = np.array([[p.r for p in row] for row in field])
mach = np.array([[p.mach for p in row] for row in field])
theta = np.array([[p.theta for p in row] for row in field])

mv.view(np.nan_to_num(x), "X")


mv.view(np.nan_to_num(r), "R")

np.nan_to_num(mach, copy=False)
mv.view(mach, "Mach")

np.nan_to_num(theta, copy=False)
mv.view(np.rad2deg(theta), "Theta")

mv.show()

Rt = np.sqrt(1 - (1/expansionRatio*np.cos(Tt)))
print(f"Throat Angle: {Tt - np.pi/2}")

thetaVx = np.cos(theta)
thetaVy = np.sin(theta)


# plt.contour(x, r, mach, levels=50, cmap='jet')
# plt.quiver(x, r, thetaVx, thetaVy, scale=25, scale_units='xy', angles='xy', headwidth=3, headlength=5, width=.002, color='black')

plt.plot(x, r, '-k', linewidth=.5) # L
plt.plot(np.transpose(x), np.transpose(r), '-k', linewidth=.5) # R
plt.plot([0, (1 - Rt)*np.tan(Tt)], [1, Rt], '-r', linewidth=2) # Throat
plt.plot(x[0,:], r[0,:], '-b', linewidth=2) # CS
plt.grid(True)
plt.show()