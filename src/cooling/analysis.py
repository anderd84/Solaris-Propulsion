from dataclasses import dataclass
import tempfile
import time
from types import NoneType
import matplotlib.pyplot as plt
import joblib
import multiprocessing as mp
from alive_progress import alive_bar, alive_it

from cooling.domain_mmap import DomainMMAP
from cooling import material, calc_cell
from general.units import Q_, unitReg
from nozzle import plots
import numpy as np

@dataclass
class DomainResistors:
    v: np.memmap | NoneType = None
    h: np.memmap | NoneType = None

def CreateDomainResistors(domain: DomainMMAP):
    workingDir = tempfile.mkdtemp()
    v = np.memmap(f'{workingDir}/vResistors.dat', dtype='float64', mode='w+', shape=(domain.vpoints - 1, domain.hpoints))
    v.fill(-1)
    h = np.memmap(f'{workingDir}/hResistors.dat', dtype='float64', mode='w+', shape=(domain.vpoints, domain.hpoints - 1))
    h.fill(-1)

    return DomainResistors(v, h)
            
def UpdateDomainResistors(domain: DomainMMAP, parallel:joblib.Parallel, resistors: DomainResistors):
    outputs = parallel(
        joblib.delayed(GetResistor)(domain, resistors, (i, j)) for i in range(domain.vpoints) for j in range(domain.hpoints)
    )

    with alive_bar(domain.vpoints*domain.hpoints, title="Precomputing resistors") as bar:
        for output in outputs:
            bar()

def GetResistor(domain: DomainMMAP, resistors, point: tuple[int, int]) -> int:
    if domain.material[point] in material.MaterialType.ADIABATIC:
        return 1
    if point[0] < domain.vpoints - 1: # vertical
        if domain.material[point[0] + 1, point[1]] in material.MaterialType.ADIABATIC:
            pass
        elif domain.material[point] in material.MaterialType.STATIC and domain.material[point[0] + 1, point[1]] in material.MaterialType.STATIC:
            pass
        elif domain.material[point] not in material.MaterialType.COOLANT or domain.material[point[0] + 1, point[1]] not in material.MaterialType.COOLANT:
            resistors.v[point] = calc_cell.CombinationResistor(domain, point, (point[0] + 1, point[1])).m_as(unitReg.hour * unitReg.degR / unitReg.BTU)

    if point[1] < domain.hpoints - 1: # horizontal
        if domain.material[point[0], point[1] + 1] in material.MaterialType.ADIABATIC:
            pass
        elif domain.material[point] in material.MaterialType.STATIC and domain.material[point[0], point[1] + 1] in material.MaterialType.STATIC:
            pass
        elif domain.material[point] not in material.MaterialType.COOLANT or domain.material[point[0], point[1] + 1] not in material.MaterialType.COOLANT:
            resistors.h[point] = calc_cell.CombinationResistor(domain, point, (point[0], point[1] + 1)).m_as(unitReg.hour * unitReg.degR / unitReg.BTU)
    return 0

def AnalyzeMC(domain: DomainMMAP, MAX_CORES: int = mp.cpu_count() - 1, tol: float = 1e-2, convPlot: bool = True, precompute: int = 0, mathmode = None):
    calcPoints = set()
    blacklist = set()
    programSolverSettings = {}
    # plt.ion()

    # fig = plots.CreateNonDimPlot()
    # imobj = domain.NodePlot(fig, "temperature", [material.DomainMaterial.CHAMBER, material.DomainMaterial.FREE])
    # fig.canvas.draw()
    # plt.pause(0.01)

    with alive_bar(domain.vpoints*domain.hpoints, title="Finding calculation points") as bar:
        for row in range(domain.vpoints):
            for col in range(domain.hpoints):
                if domain.material[row, col] not in material.MaterialType.STATIC_TEMP:
                    calcPoints.add((row, col))
                if domain.material[row, col] == material.DomainMaterial.COOLANT_BULK and domain.border[row, col]:
                    blacklist.add(tuple(domain.previousFlow[row, col]))
                bar()
    
    for pair in blacklist:
        calcPoints.remove(pair)

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
        if precompute:
            print("Startng precompute")
            resistorMap = CreateDomainResistors(domain)
            # UpdateDomainResistors(domain, parallel, res)
            programSolverSettings['resistors'] = resistorMap
        if mathmode is not None:
            programSolverSettings['mathmode'] = mathmode

        print(f" Starting solve")
        print(f" Settings: ")
        print(f"   Max Cores: {MAX_CORES}")
        print(f"   Precompute: {precompute}")
        print(f"   Tolerance: {tol}")
        print(f"   Mathmode: {mathmode}")
        while diff > tol:
            tic = time.time()
            i += 1
            diff = 0
            maxT = 0
            changes = []

            # if precompute:# and i % 5 == 0:
            #     print("Updating precompute")
            #     UpdateDomainResistors(domain, parallel, res)
            #     # programSolverSettings['resistors'] = res

            
            with alive_bar(len(calcPoints), title=f"Analyzing iteration {i}") as bar:
                outputs = parallel(
                    joblib.delayed(calc_cell.CalculateCell)(domain, row, col, **programSolverSettings) for row, col in calcPoints
                )

                for output in outputs:
                    for changeOrder in output:
                        changes.append(changeOrder)
                        row, col = changeOrder.row, changeOrder.col
                        if changeOrder.temperature is not None:
                            newTemp = changeOrder.temperature
                            newDiff = abs(domain.temperature[row, col].magnitude - newTemp.to(domain.units["temperature"]).magnitude) / domain.temperature[row, col].magnitude
                            diff = max(diff, newDiff)
                            maxT = max(maxT, newTemp.magnitude)
                    bar()

            if precompute and i % precompute == 0:
                print("reseting precomputed resistors")
                resistorMap.v.fill(-1)
                resistorMap.h.fill(-1)

            for changeOrder in alive_it(changes):
                row, col = changeOrder.row, changeOrder.col
                # if precompute and i % precompute == 0:
                    # resistorMap.v[[min(row, domain.vpoints-2) ,row-1], col] = -1
                    # resistorMap.h[row, [min(col, domain.hpoints - 2), col-1]] = -1

                if changeOrder.temperature is not None:
                    currentTemp = domain.temperature[row, col]
                    newTemp = changeOrder.temperature
                    diff_ = newTemp.m - currentTemp.m
                    # maxChange = currentTemp.m / 100
                    maxChange = 2 * np.sign(diff_)
                    if False and (domain.material[row, col] in material.MaterialType.COOLANT or domain.border[row,col]):
                        domain.setMEM(row, col, 'temperature', Q_(max(min(abs(diff_), maxChange), -abs(diff_)),currentTemp.units) + currentTemp)
                        # if abs(diff_) > 100:
                        #     domain.setMem(row, col, 'temperature', currentTemp - diff)
                    else:
                        domain.setMEM(row, col, 'temperature', changeOrder.temperature)
                if changeOrder.pressure is not None:
                    domain.setMEM(row, col, 'pressure', changeOrder.pressure)

            toc = time.time()
            print(f"Total computation time: {toc - tic}")

            print(f"Max diff: {diff*100}%")
            print(f"Max temp: {maxT}R")
            diffArr.append(diff*100)

            # imobj = domain.updateNodePlot(imobj, "temperature", [material.DomainMaterial.CHAMBER, material.DomainMaterial.FREE])
            # fig.canvas.draw()
            # fig.canvas.flush_events()
            # plt.pause(0.01)

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

            if i % 50 == 0:
                print("saving progress")
                mesh = domain.toDomain()
                mesh.DumpFile("save")
        
            # if maxT > 1000:
            #     print("max temp too high, stopping")
            #     mesh = domain.toDomain()
            #     mesh.DumpFile("save")
            #     break

            if i == 500:
                break

    # if precompute:
    #     print("plotting relations")
    #     for row in range(domain.vpoints):
    #         for col in range(domain.hpoints):
    #             if row < domain.vpoints - 1 and resistorMap.v[row, col] != -1:
    #                 c = 'g' if resistorMap.v[row, col] > 0 else 'r'
    #                 c = 'k' if np.isclose(resistorMap.v[row, col], 0) else c
    #                 plt.plot([domain.x[row, col], domain.x[row + 1, col]], [domain.r[row, col], domain.r[row + 1, col]], c, linewidth=1)
    #             if col < domain.hpoints - 1 and resistorMap.h[row, col] != -1:
    #                 c = 'g' if resistorMap.h[row, col] > 0 else 'r'
    #                 c = 'k' if np.isclose(resistorMap.h[row, col], 0) else c
    #                 plt.plot([domain.x[row, col], domain.x[row, col + 1]], [domain.r[row, col], domain.r[row, col + 1]], c, linewidth=1)
    #     np.product = np.prod
    #     import matrix_viewer as mv
    #     mv.view(resistorMap.v)
    #     mv.view(resistorMap.h)
    #     mv.show_with_pyplot()

    print("saving progress")
    mesh = domain.toDomain()
    mesh.DumpFile("save" + str(mathmode))
    print("-------------------------")
    print(f" Done! in {i} iterations ")
    print("-------------------------")