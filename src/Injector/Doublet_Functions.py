import numpy as np
def spike_contour(Points):
    # -------------- Code to Make my Shitty version of the Spike COntour -------------- #
    #Constants for Spike contour
    r1 = 3.50
    r2 = 2.50
    angle1 = 41.81
    angle2 = 83.62
    # Initial positions
    startX, startY = 0, 1
    startX1, startY1 = 3.5, 1  # Adjusted based on given code
    BaseX = np.linspace(startX, startX1, Points)
    BaseY = np.linspace(startY, startY1, Points)
    # Calculating end points and start points for the arcs
    endX1 = startX1 + r1 * np.sin(np.radians(angle1))
    endY1 = startY1 + r1 * (1 - np.cos(np.radians(angle1)))
    startX2, startY2 = endX1, endY1
    endX2 = startX2 + r2 * (np.cos(np.radians(90 - angle2 / 2)) - np.cos(np.radians(90 + angle2 / 2)))
    endY2 = startY2 + r2 * (np.sin(np.radians(90 - angle2 / 2)) - np.sin(np.radians(90 + angle2 / 2)))
    startX3, startY3 = endX2, endY2
    endX3 = startX3 + r1 * (np.cos(3*np.pi/2) - np.cos(np.radians(270 - angle1)))
    endY3 = startY3 + r1 * (np.sin(3*np.pi/2) - np.sin(np.radians(270 - angle1)))
    # Defining theta ranges for the arcs
    thetaRange1 = np.linspace(0, np.radians(angle1), Points)
    thetaRange2 = np.linspace(np.radians(90 + angle2 / 2), np.radians(90 - angle2 / 2), Points)
    thetaRange3 = np.linspace(np.radians(angle1), 0, Points)
    # Defining arcs
    arc1_x = startX1 + r1 * np.sin(thetaRange1)
    arc1_y = startY1 + r1 * (1 - np.cos(thetaRange1))
    arc2_x = startX2 + r2 * (np.cos(np.radians(90 - angle2 / 2)) + np.cos(thetaRange2))
    arc2_y = startY2 + r2 * (-np.sin(np.radians(90 - angle2 / 2)) + np.sin(thetaRange2))
    arc3_x = endX3 + r1 * np.cos(3 * np.pi / 2 - thetaRange3)
    arc3_y = endY3 + r1 * (1 + np.sin(3 * np.pi / 2 - thetaRange3))
    # Combining all segments (for both plotting use and finding peaks)
    x_profile = np.concatenate([BaseX, arc1_x, arc2_x, arc3_x])
    y_profile = np.concatenate([BaseY, arc1_y, arc2_y, arc3_y]) + 0.5
    return x_profile,y_profile