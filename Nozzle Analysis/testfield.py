import rao
import numpy as np
import matplotlib.pyplot as plt
import matrix_viewer as mv

gamma = 1.23
Me = 2.4
Te = np.deg2rad(-8.25)
Mt = 1.15
PbPc = 0

areaRatio, Cf, lengthRatio, aplot, rplot, lplot, tplot, mplot = rao.InputDataGenerate(Me, Te, gamma, .00001, PbPc)

print(areaRatio, Cf, lengthRatio)
# input("Press Enter to continue...")


Tt = rao.CalculateThroatAngle(Me, Te, Mt, gamma)
# length = 1.9554
# expansionRatio = 5

length = lengthRatio
expansionRatio = areaRatio

controlSurface: np.ndarray[rao.CharacteristicPoint] = rao.GetControlSurfaceProperties(Me, Te, length, gamma, 50)
expansionFan: np.ndarray[rao.CharacteristicPoint] = rao.GenerateExpansionFan(Me, Mt, Tt, gamma , 50)

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


Rt = np.sqrt(1 - (1/expansionRatio*np.cos(Tt)))
print(f"Throat Angle: {Tt - np.pi/2}")

cont = rao.CalculateContour(field, Rt, Tt)
cx = [p.x for p in cont]
cy = [p.r for p in cont]

field = rao.PruneField(field)

field = rao.PruneUnderContour(field, cont)


x = np.array([[p.x for p in row] for row in field])
r = np.array([[p.r for p in row] for row in field])
mach = np.array([[p.mach for p in row] for row in field])
theta = np.array([[p.theta for p in row] for row in field])
alpha = np.array([[p.alpha for p in row] for row in field])

mv.view(np.nan_to_num(x), "X")

mv.view(np.nan_to_num(r), "R")

mv.view(np.nan_to_num(alpha), "alpha")

np.nan_to_num(mach, copy=False)
mv.view(mach, "Mach")

np.nan_to_num(theta, copy=False)
mv.view(np.rad2deg(theta), "Theta")

mv.view(np.transpose(np.array([cx, cy])), "Contour")

mv.show()

csarrows = 15
fanarrows = 5

qx = x[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows]
qy = r[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows]


thetaVx = np.cos(theta[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows])
thetaVy = np.sin(theta[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows])

plt.contourf(x, r, mach, levels=100, cmap='jet')
plt.legend(['Mach'])
plt.quiver(qx, qy, thetaVx, thetaVy, scale=25, scale_units='xy', angles='xy', headwidth=3, headlength=5, width=.002, color='black')


plt.plot(x, r, '-k', linewidth=.5) # L
plt.plot(np.transpose(x), np.transpose(r), '-k', linewidth=.5) # R
plt.plot([0, (1 - Rt)*np.tan(Tt)], [1, Rt], '-r', linewidth=2) # Throat
plt.plot(x[0,:], r[0,:], '-b', linewidth=2) # CS
plt.grid('on', linestyle='--')
plt.gca().set_xlabel('X/Re')
plt.gca().set_ylabel('R/Re')

cx.append((1 - Rt)*np.tan(Tt))
cx.append(x[0,-1])
cy.append(0)
cy.append(0)

plt.fill(cx, cy, 'k')

plt.show()

# pt1 = field[0, -2]
# pt2 = field[1, -1]
# L = field[0, -1]

# cont = [L]

# m = 0
# n = 0

# rows = field.shape[0]
# cols = field.shape[1]

# def eqn(pt1: rao.CharacteristicPoint, pt2: rao.CharacteristicPoint, L: rao.CharacteristicPoint, ax: rao.CharacteristicPoint):
#     return (L.r - pt1.r + (pt1.x - ax)*(pt1.r-pt2.r)/(pt1.x - pt2.x))/(L.x - ax) - np.tan((pt1.x-ax)/(pt1.x-pt2.x)*pt2.theta+(ax-pt2.x)/(pt1.x-pt2.x)*pt1.theta)

# def calcNext(pt1, pt2, L):
#     ax2 = fsolve(lambda ax: eqn(pt1, pt2, L, ax), pt2.x)[0]
#     ar2 = pt1.r - (pt1.x - ax2)*(pt1.r - pt2.r)/(pt1.x - pt2.x)
#     print(ax2, ar2)
#     return rao.CharacteristicPoint(ax2, ar2, 0, 0)

# plt.ion()
# points, = plt.plot([pt1.x, pt2.x], [pt1.r, pt2.r], '.k', markersize=10) # Points
# xl = np.array([p.x for p in cont])
# rl = np.array([p.r for p in cont])

# line0, = plt.plot(xl, rl, '-k', linewidth=2) # Contour

# while cont[-1].x not in x[:, 1]:
#     pt1 = field[0 + m, -2 + n]
#     pt2 = field[1 + m, -1 + n]
#     try:
#         Lnew = calcNext(pt1, pt2, cont[-1]) # 1
#     except:
#         print(f"Failed at {(m, n)}")
#         break
#     if Lnew.x > max(pt2.x, pt1.x) or Lnew.x < min(pt2.x, pt1.x):
#         print(f"pt1.x = {pt1.x}, pt2.x = {pt2.x}, Lnew.x = {Lnew.x}")
#         m += 1
#         n += 1
        
#     else:
#         cont.append(Lnew)
#         n -= 1
#     if m >= rows - 1 or n <= -(cols - 1):
#         print(f"pt1.x = {pt1.x}, pt2.x = {pt2.x}, Lnew.x = {Lnew.x}")
#         break

#     points.set_xdata([pt1.x, pt2.x])
#     points.set_ydata([pt1.r, pt2.r])

#     xl = np.array([p.x for p in cont])
#     rl = np.array([p.r for p in cont])
#     line0.set_xdata(xl)
#     line0.set_ydata(rl)

#     plt.draw()
#     plt.pause(.001)
#     # plt.waitforbuttonpress()

# xl = np.append(xl, (1 - Rt)*np.tan(Tt))
# rl = np.append(rl, Rt)
# plt.plot(xl, rl, '-k', linewidth=2) # Contour
