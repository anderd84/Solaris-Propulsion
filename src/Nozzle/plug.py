import numpy as np
import matplotlib.pyplot as plt

from fluids.gas import Gas






def CreateRaoPlug(exhaustGas: Gas, chamberPressure: float, designAmbient: float, lipRadius: float, maxLengthConsidered: float):
    exhaustGas.P0 = chamberPressure
    