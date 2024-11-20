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

def AnalyzeMC(domain: domain.DomainMMAP, fig, plugC, cowlC, MAX_CORES: int = mp.cpu_count() - 1, tol: float = 1e-2, convPlot: bool = True):
    # plt.show()

    # for col in range(domain.hpoints):
    #     for row in range(domain.vpoints):
    #         outs = analysisCoolingRef.Cell(domain, row, col)
    #         domain.setMEM(row, col, 'temperature', outs[0])
    #         domain.setMEM(row, col, 'pressure', outs[1])
            # if domain.material[row, col] not in material.MaterialType.STATIC_TEMP:
            #     fig.axes[0].clear()
            #     plots.PlotPlug(fig, plugC)
            #     plots.PlotPlug(fig, cowlC)
            #     domain.toDomain().ShowStatePlot(fig)
            #     fig.axes[0].plot([domain.x[row,col]], [domain.r[row,col]], 'wx')
            #     print(outs)

            #     fig.canvas.draw()
            #     fig.canvas.flush_events()
            #     plt.waitforbuttonpress()
            #     input("Press Enter to continue...")
    if convPlot:
        plt.ion()
        convergePlot, ax = plt.subplots()
        ax.set_title("Convergence")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Max % Difference")
        ax.grid(True)
        convergeP, = ax.plot([1], [1], 'r-')

    diffArr = []

    with joblib.Parallel(n_jobs=MAX_CORES, return_as='generator') as parallel:
        diff = tol + 1
        numRows = domain.vpoints
        i = 0
        while diff > tol:
            i += 1
            diff = 0
            with alive_bar(domain.vpoints, title=f"Analyzing iteration {i}") as bar:
                outputs = parallel(
                    joblib.delayed(CalcRow)(domain, row) for row in range(numRows)
                )

                for output in outputs:
                    for row, col, tp in output:
                        diff = max(diff, abs(domain.temperature[row, col].magnitude - tp[0].to(domain.units["temperature"]).magnitude)/domain.temperature[row, col].magnitude)
                        domain.setMEM(row, col, 'temperature', tp[0])
                        domain.setMEM(row, col, 'pressure', tp[1])
                    bar()

            print(f"Max diff: {diff*100}%")
            diffArr.append(diff*100)
            if convPlot:
                del convergeP
                convergeP = ax.plot(range(i), diffArr, 'r-')
                ax.autoscale(axis='y')
                
                convergePlot.canvas.draw()
                convergePlot.canvas.flush_events()

            if i % 10 == 0:
                print("saving progress")
                mesh = domain.toDomain()
                mesh.DumpFile("save.msh")





def CalcRow(domain: domain.DomainMMAP, row: int):
    res = []
    for col in range(domain.hpoints):
        res.append((row, col, analysisCoolingRef.Cell(domain, row, col)))
    return res