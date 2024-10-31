import numpy as np
from icecream import ic
import matrix_viewer as mv

thetaL = np.deg2rad(15)
thetal = np.deg2rad(90 - 10 - abs(68.5157))

ic(np.rad2deg(thetal))

sL = np.sin(thetaL)
cL = np.cos(thetaL)
sl = np.sin(thetal)
cl = np.cos(thetal)
cotl = 1/np.tan(thetal)
cscl = 1/np.sin(thetal)

o = .025
Re = 3
Rinner = 3.5
Rmax = 4
xc = -.5456
xm = .15

Amat = np.array([[0, 0, 1, 0, sL, 0, 0, 0, 0], 
                [0, 0, 0, 1, -cL, 0, 0, 0 ,0],
                [1, 1, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 1, 1, 0, 0],
                [0, -cl, 1, 0, -cl, 0, 0, -sl, 0],
                [1, sl, 0, -1, sl, 0, 0, -cl, 0],
                [0, 0, 0, 0, 0, 0, -sL, 0, cL],
                [0, 0, 0, 0, 0, 1, -cL, 0, -sL],
                [0, 0, cotl, 1, -cscl, 0, 0, 0, 0]])

bmat = np.array([0, Re, Rinner, Rmax, xc, 0, xm, Re, Re-o*cscl])

ic(np.linalg.solve(Amat, bmat))