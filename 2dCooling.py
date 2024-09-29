from matplotlib import pyplot as plt, patches
from Doublet_Functions import spike_contour
from InjectorCad import injector_cad_write
from Drill import drill_approximation
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from fluid import Q_, ureg, CD_drill, Pressure_Drop_Fuel, Pressure_Drop_Lox,  pm
from Doublet import OX_CORE, FUEL_CORE
import fluid
import numpy as np
import os

#First step always is to update doublet.py file and run before hand to grab all mdot and density values at injector side





