import numpy as np
import matplotlib.pyplot as plt
from icecream import ic

from fluids import gas
import general.design as DESIGN
from general.units import Q_, unitReg
from nozzle import angelino

exhaust = DESIGN.exhaustGas
astar = DESIGN.chokeArea
ic(astar)
angelino.CalculateInternalExternalNozzle(2,2.916,.15,exhaust, 100)