import matplotlib.pyplot as plt
import joblib
import multiprocessing as mp
from alive_progress import alive_bar

from cooling.domain import DomainMC, DomainMMAP, SparseDomain
from cooling import material, calc_cell

def AnalyzeMC(domain: DomainMMAP, MAX_CORES: int = mp.cpu_count() - 1, tol: float = 1e-2, convPlot: bool = True):
    calcPoints = []
    with alive_bar(domain.vpoints*domain.hpoints, title="Finding calculation points") as bar:
        for row in range(domain.vpoints):
            for col in range(domain.hpoints):
                if domain.material[row, col] not in material.MaterialType.STATIC_TEMP:
                    calcPoints.append((row, col))
                bar()

    if convPlot:
        plt.ion()
        convergePlot, ax = plt.subplots()
        ax.set_title("Convergence")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Max % Difference")
        ax.grid(True)

    diffArr = []

    with joblib.Parallel(n_jobs=MAX_CORES, return_as='generator') as parallel:
        diff = tol + 1
        numRows = domain.vpoints
        i = 0
        while diff > tol:
            i += 1
            diff = 0
            with alive_bar(len(calcPoints), title=f"Analyzing iteration {i}") as bar:
                outputs = parallel(
                    joblib.delayed(CalcCell)(domain, row, col) for row, col in calcPoints
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
                # del convergeP
                ax.clear()
                ax.grid(True)
                iarr = max(0, i - 15)
                ax.plot(range(iarr, i), diffArr[iarr:], 'k-')
                # convergeP = ax.plot(range(i), diffArr, 'r-')
                ax.autoscale(axis='y')
                
                convergePlot.canvas.draw()
                convergePlot.canvas.flush_events()

            if i % 10 == 0:
                print("saving progress")
                mesh = domain.toDomain()
                mesh.DumpFile("save")
            
            if i == 2:
                break

def AnalyzeMCSparse(domain: DomainMC, MAX_CORES: int = mp.cpu_count() - 1, tol: float = 1e-2, convPlot: bool = True):
    calcPoints = []
    with alive_bar(domain.vpoints*domain.hpoints, title="Finding calculation points") as bar:
        for row in range(domain.vpoints):
            for col in range(domain.hpoints):
                if domain.array[row, col].material not in material.MaterialType.STATIC_TEMP:
                    calcPoints.append((row, col))
                bar()

    if convPlot:
        plt.ion()
        convergePlot, ax = plt.subplots()
        ax.set_title("Convergence")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Max % Difference")
        ax.grid(True)

    diffArr = []

    with joblib.Parallel(n_jobs=MAX_CORES, return_as='generator') as parallel:
        diff = tol + 1
        i = 0
        while diff > tol:
            i += 1
            diff = 0
            with alive_bar(len(calcPoints), title=f"Analyzing iteration {i}") as bar:
                outputs = parallel(
                    joblib.delayed(CalcCell)(SparseDomain(domain, (row, col)), row, col) for row, col in calcPoints
                )

                for output in outputs:
                    for row, col, tp in output:
                        diff = max(diff, abs(domain.array[row, col].temperature.magnitude - tp[0].magnitude)/domain.array[row, col].temperature.magnitude)
                        domain.array[row,col].temperature = tp[0]
                        domain.array[row,col].pressure = tp[1]
                    bar()

            print(f"Max diff: {diff*100}%")
            diffArr.append(diff*100)
            if convPlot:
                # del convergeP
                ax.clear()
                ax.grid(True)
                iarr = max(0, i - 15)
                ax.plot(range(iarr, i), diffArr[iarr:], 'k-')
                # convergeP = ax.plot(range(i), diffArr, 'r-')
                ax.autoscale(axis='y')
                
                convergePlot.canvas.draw()
                convergePlot.canvas.flush_events()

            if i % 10 == 0:
                print("saving progress")
                domain.DumpFile("save")

            if i == 2:
                break

def CalcCell(domain: DomainMMAP | SparseDomain, row: int, col: int):
    return [(row, col, calc_cell.Cell(domain, row, col))]