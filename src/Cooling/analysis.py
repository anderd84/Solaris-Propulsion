# ok this is gonna be complicated

from Nozzle import plots
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
import joblib
import multiprocessing as mp
from alive_progress import alive_bar, alive_it
from General.units import Q_, unitReg

from Cooling import domain
from Cooling import material
from Cooling import analysisCoolingRef

def AnalyzeMC(domain: domain.DomainMMAP, MAX_CORES: int = mp.cpu_count() - 1, tol: float = 1e-2):
    with joblib.Parallel(n_jobs=MAX_CORES, return_as='generator') as parallel:
        diff = tol + 1
        numRows = domain.vpoints
        for i in range(1):
            with alive_bar(domain.vpoints, title=f"Analyzing iteration {i}") as bar:
                outputs = parallel(
                    joblib.delayed(CalcRow)(domain, row) for row in range(numRows)
                )

                for output in outputs:
                    for row, col, tp in output:
                        domain.setMEM(row, col, 'temperature', tp[0])
                        domain.setMEM(row, col, 'pressure', tp[1])
                    bar()





def CalcRow(domain: domain.DomainMMAP, row: int):
    res = []
    for col in range(domain.hpoints):
        res.append((row, col, analysisCoolingRef.Cell(domain, row, col)))
    return res