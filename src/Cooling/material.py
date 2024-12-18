
from matplotlib import pyplot as plt, patches
import numpy as np
from icecream import ic
from Nozzle.nozzle import ContourPoint
from enum import Enum

class DomainMaterial(Enum):
    FREE = 0 
    COWL = 1
    COOLANT = 2
    CHAMBER = 3
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