from matplotlib import pyplot as plt, patches
import numpy as np
import os
from scipy.optimize import fsolve
from icecream import ic

from fluids.fluid import CD_drill, Pressure_Drop_Fuel, \
                         Pressure_Drop_Lox,  pm, get_fluid_properties, CP
from General.units import Q_, unitReg
from Cooling.cooling2d import combustion_convection
from Cooling.domain import DomainPoint

def channel_area_sizing(point: DomainPoint):
    # Chooses a channel area based on the temperature of the inner nozzle wall
    # T_wi = inner wall temperature
    # A_c = current cross-sectional area of cooling channel 

    # Target a certain inner wall temperature
    # If above, decrease area to increase h_cool
    # Otherwise, increase area to reduce pressure losses

    # Lame implementation
    T_wi = point.temperature
    T_wi = T_wi.to(T_wi, unitReg.degR)
    if T_wi >= 3400:
        point.area -= 0.0005
    if T_wi < 3000:
        point.area += 0.0005
