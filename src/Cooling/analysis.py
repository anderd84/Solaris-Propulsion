# ok this is gonna be complicated

import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
import joblib
import multiprocessing as mp

from Cooling import domain
from Cooling import material

def AnalyzeMC(domain: domain.DomainMMAP, MAXCORES: int = mp.cpu_count() - 1, tol: float = 1e-2):
    with joblib.Parallel(n_jobs=MAXCORES) as parallel:
        diff = tol + 1
        numRows = domain.vpoints
        while diff > tol:
            outputs = parallel(
                joblib.delayed(CalcRow)(domain, row) for row in range(numRows)
            )

def CalcRow(domain: domain.DomainMMAP, row: int):
    pass