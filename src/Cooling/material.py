
from matplotlib import pyplot as plt, patches
import numpy as np
from icecream import ic
from Nozzle.nozzle import ContourPoint
from enum import Enum

class DomainMaterial(Enum):
    FREE = 0 
    COWL = 1
    COOLANT = 2
    CHAMBER = 3 #*Gas inside the chamber
    PLUG = 4

def isIntersect(Point, contour: np.ndarray[ContourPoint], domainSize: tuple[int, int]):
    if contour.size == 0:
        return False
    OnTheLine = 0
    A = [Point.x, Point.r]   #Point
    B = [Point.x, Point.r + domainSize[1]]   #Making an "infinitely" long vertical line at the point
    Check = 0
    for i in range(len(contour) - 1):
        C = [contour[i].x,contour[i].r]
        D = [contour[i+1].x,contour[i+1].r]
        Check = Check + (ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D))
        if is_point_on_segment(Point, C, D):
            OnTheLine = 1
            break
    C = [contour[-1].x,contour[-1].r]
    D = [contour[0].x,contour[0].r]
    Check = Check + (ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D))    

    InsideOutside = np.mod(Check,2) == 1 or OnTheLine == 1
    return InsideOutside #Return true if point is inside the polygon

def is_point_on_segment(Point, C, D):
    Point = [Point.x, Point.r]
    # Check collinearity using cross product
    collinear = (D[1] - C[1]) * (Point[0] - C[0]) == (Point[1] - C[1]) * (D[0] - C[0])
    # Check if point P is within the bounding box of C and D
    within_bounds = min(C[0], D[0]) <= Point[0] <= max(C[0], D[0]) and min(C[1], D[1]) <= Point[1] <= max(C[1], D[1])
    return collinear and within_bounds

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def intersectPolyAt(polygon, point1, point2):
    Sx, Sy = point1
    Tx, Ty = point2
    for j in range(len(polygon) - 1):
        a, b, c, d = polygon[j].x, polygon[j].r, polygon[j+1].x, polygon[j+1].r
        Tx = Tx if Tx != Sx else Tx - 1e-6
        c = c if c != a else c - 1e-6
        mL = (Ty - Sy)/(Tx - Sx)
        mC = (d - b)/(c - a)
        Amat = np.array([[mL, -1], [mC, -1]])
        bmat = np.array([[mL*Sx - Sy], [mC*a - b]])
        X = np.linalg.solve(Amat, bmat)
        Bx = X[0,0]
        By = X[1,0]
        if min(Sx,Tx) <= Bx <= max(Sx,Tx) and min(Sy,Ty) <= By <= max(Sy,Ty) and min(a,c) <= Bx <= max(a,c) and min(b,d) <= By <= max(b,d):
            return (Bx, By), (j, j+1)
    return None, None

def main():
    Points = 1000

    contour_x = np.concatenate([np.linspace(0, 5, Points),np.linspace(5, 5, Points),np.linspace(5, 0, Points),np.linspace(0, 0, Points)])
    contour_r = np.concatenate([np.linspace(0, 0, Points),np.linspace(0, 5, Points),np.linspace(5, 5, Points),np.linspace(5, 0, Points)])

    contour = [ContourPoint(row[0],row[1]) for row in np.transpose([contour_x,contour_r])]
    Point = ContourPoint(2.5,2.5)

    Check = isIntersect(Point, contour, np.array([10000,10000]))
    ic(Check)







if __name__ == "__main__":
    main()