import matplotlib.pyplot as plt
import numpy as np
from Nozzle.nozzle import ContourPoint
from Nozzle.rao import CharacteristicPoint
from scipy.optimize import fsolve
import copy

def CreateNonDimPlot() -> plt.Figure:
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel('X/Re')
    ax.set_ylabel('R/Re')
    return fig

def PlotField(fig: plt.Figure, field: np.ndarray, scale = 1, csarrows: int = 15, fanarrows: int = 10) -> plt.Figure:
    x = np.array([[p.x*scale for p in row] for row in field])
    r = np.array([[p.r*scale for p in row] for row in field])
    mach = np.array([[p.mach for p in row] for row in field])
    theta = np.array([[p.theta for p in row] for row in field])
    alpha = np.array([[p.alpha for p in row] for row in field])

    qx = x[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows]
    qy = r[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows]


    thetaVx = np.cos(theta[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows])
    thetaVy = np.sin(theta[::field.shape[0]//fanarrows, ::field.shape[1]//csarrows])

    ax = fig.axes[0]

    # machContours = ax.contourf(x, r, mach, levels=100, cmap='jet')
    # fig.colorbar(machContours, orientation='vertical')
    ax.quiver(qx, qy, thetaVx, thetaVy, scale=25, scale_units='xy', angles='xy', headwidth=3, headlength=5, width=.002, color='black')

    ax.plot(x, r, '-k', linewidth=.5) # L
    ax.plot(np.transpose(x), np.transpose(r), '-k', linewidth=.5) # R
    ax.plot(x[0,:], r[0,:], '-b', linewidth=2) # CS
    ax.grid('on', linestyle='--')

    return fig

def PlotContour(fig: plt.Figure, contour: np.ndarray[ContourPoint | CharacteristicPoint], Rt, Tt, lipRadius = 1) -> plt.Figure:
    ax = fig.axes[0]

    cx = [p.x for p in contour]
    cy = [p.r for p in contour]

    if lipRadius == 1:
        cx.append((lipRadius - Rt)*np.tan(Tt))
        cy.append(0)

        cx.append(cx[0])
        cy.append(0)

        cx.append(cx[0])
        cy.append(cy[0])

    ax.plot([0, (lipRadius - Rt)*np.tan(Tt)], [lipRadius, Rt], '-r', linewidth=2) # Throat
    ax.plot(cx, cy, '-k', linewidth=2) # Contour
    # ax.fill(cx, cy, 'k')
    
    return fig

def PlotPlug(fig: plt.Figure, plug: np.ndarray[ContourPoint]) -> plt.Figure:
    ax = fig.axes[0]

    x = [p.x for p in plug]
    r = [p.r for p in plug]

    ax.plot(x, r, '-k', linewidth=2)

    return fig

def show3d(contour: np.ndarray):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    x = np.array([p.x for p in contour])
    r = np.array([p.r for p in contour])

    thetas = np.linspace(0, 2*np.pi, 101)

    y = np.outer(np.cos(thetas), r)
    z = np.outer(np.sin(thetas), r)
    x = np.outer(np.ones(len(thetas)), x)

    ax.plot_surface(x, y, z, rstride=1, cstride=1, color='b', alpha=0.5)
    ax.set_xlabel('X')
    ax.set_ylabel('R')
    ax.set_zlabel('Theta')

    plt.show()

def WriteContourTXT(contour: np.ndarray[CharacteristicPoint], filename: str):
    with open(filename, 'w') as f:
        for point in contour:
            f.write(f'{point.x},{point.r},0\n')

def LiveContour(field: np.ndarray[CharacteristicPoint], Rt: float, Tt: float, plot: plt.Figure) -> np.ndarray[CharacteristicPoint]:
    ax = plot.axes[0]

    plt.ion()

    pt1 = field[0, -2]
    pt2 = field[1, -1]
    L = field[0, -1]

    cont = [L]

    m = 0
    n = 0

    rows, cols = field.shape

    pts, = ax.plot([pt1.x, pt2.x], [pt1.r, pt2.r], '.w', markersize=10)
    contLine, = ax.plot([L.x], [L.r], '-w', linewidth=2)

    Lold = L.clone()

    def eqn(pt1: CharacteristicPoint, pt2: CharacteristicPoint, L: CharacteristicPoint, ax: CharacteristicPoint) -> float:
        return (L.r - pt1.r + (pt1.x - ax)*(pt1.r-pt2.r)/(pt1.x - pt2.x))/(L.x - ax) - np.tan((pt1.x-ax)/(pt1.x-pt2.x)*pt2.theta+(ax-pt2.x)/(pt1.x-pt2.x)*pt1.theta)

    def calcNext(pt1: CharacteristicPoint, pt2: CharacteristicPoint, L: CharacteristicPoint) -> CharacteristicPoint:
        ax2 = fsolve(lambda ax: eqn(pt1, pt2, L, ax), pt2.x)[0]
        ar2 = pt1.r - (pt1.x - ax2)*(pt1.r - pt2.r)/(pt1.x - pt2.x)
        # print(ax2, ar2)
        return CharacteristicPoint(ax2, ar2, np.arctan((ar2 - L.r)/(ax2 - L.x)), 0)

    while m < rows - 2 and n > -(cols - 2):
        pt1 = field[0 + m, -2 + n]
        pt2 = field[1 + m, -1 + n]

        try:
            Lnew = calcNext(pt1, pt2, cont[-1]) # 1
        except:
            print(f"Failed at {(m, n)}")
            break

        pts.set_xdata([pt1.x, pt2.x, Lnew.x])
        pts.set_ydata([pt1.r, pt2.r, Lnew.r])

        contLine.set_xdata([p.x for p in cont])
        contLine.set_ydata([p.r for p in cont])

        plot.canvas.draw()
        plot.canvas.flush_events()
        print(f"m: {m}, n: {n}")
        plt.waitforbuttonpress()

        if Lnew.r < 0:
            print("Negative R")
            cont.append(Lold)
            n -= 2
            m -= 1
        else:
            if Lnew.x > max(pt2.x, pt1.x) or Lnew.x < min(pt2.x, pt1.x):
                m += 1
                n += 1
                Lold = copy.copy(Lnew)
            else:
                cont.append(Lnew)
                n -= 1


    cont.append(CharacteristicPoint((1 - Rt)*np.tan(Tt), Rt, 0, 0))

    return np.array(cont)