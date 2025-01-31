from types import NoneType
import cooling
import numpy as np
import pint
from dataclasses import dataclass

from cooling import cooling2d, domain, domain_mmap
from cooling.material import DomainMaterial, MaterialType
from fluids import fluid
import general.design as DESIGN
from general.units import Direction, unitReg, Q_

class CellUpdates:
    row: int
    col: int
    temperature: pint.Quantity | NoneType
    pressure: pint.Quantity | NoneType

    def __init__(self, row, col, **kwargs):
        self.row = row
        self.col = col
        self.temperature = kwargs.get("temperature", None)
        self.pressure = kwargs.get("pressure", None)

@dataclass
class ThermalResistor:
    R: pint.Quantity = Q_(-1, unitReg.hour * unitReg.degR / unitReg.BTU)
    T: pint.Quantity = Q_(-1, unitReg.degR)

def getConductivity(domain: domain_mmap.DomainMMAP, row: int, col: int) -> pint.Quantity:
    if domain.material[row,col] in MaterialType.WALL:
        return cooling2d.conduction_grcop(domain.temperature[row,col].to(unitReg.degR))
    else:
        return cooling2d.conduction_rp1(domain.temperature[row,col].to(unitReg.degR))

def ConductionResistor(domain: domain_mmap.DomainMMAP, sink: tuple[int, int], source: tuple[int, int]) -> pint.Quantity: # sink is to, source is from
    L = Q_(domain.xstep, unitReg.inch).to(unitReg.foot)
    k = getConductivity(domain, sink[0], sink[1])
    if source[0] == sink[0]: # same row, horizontal conduction
        # conductionArea = (2 * np.pi * Q_(domain.r[sink], unitReg.inch).to(unitReg.foot) ) * Q_(domain.rstep, unitReg.inch).to(unitReg.foot)
        ro = Q_(domain.r[sink] + domain.rstep/2, unitReg.inch).to(unitReg.foot)
        ri = Q_(domain.r[sink] - domain.rstep/2, unitReg.inch).to(unitReg.foot)
        conductionArea = np.pi*np.abs(ro**2 - ri**2)
        return L / (conductionArea * k)
    elif source[1] == sink[1]: # same column, vertical conduction
        ri = Q_(min(domain.r[source], domain.r[sink]), unitReg.inch).to(unitReg.foot)
        ro = Q_(max(domain.r[source], domain.r[sink]), unitReg.inch).to(unitReg.foot)
        return np.log(ro / ri) / (2 * np.pi * L * k)

def ConductionHalfResistor(domain: domain_mmap.DomainMMAP, sink: tuple[int, int], source: tuple[int, int], sinkSide: bool = True) -> pint.Quantity:
    L = Q_(domain.xstep, unitReg.inch).to(unitReg.foot)
    k = getConductivity(domain, sink[0], sink[1]) if sinkSide else getConductivity(domain, source[0], source[1])
    if source[0] == sink[0]: # same row, horizontal conduction
        # conductionArea = (2 * np.pi * Q_(domain.r[sink], unitReg.inch).to(unitReg.foot)) * Q_(domain.rstep, unitReg.inch).to(unitReg.foot)
        ro = Q_(domain.r[sink] + domain.rstep/2, unitReg.inch).to(unitReg.foot)
        ri = Q_(domain.r[sink] - domain.rstep/2, unitReg.inch).to(unitReg.foot)
        conductionArea = np.pi*np.abs(ro**2 - ri**2)
        return (L / 2) / (conductionArea * k)
    elif source[1] == sink[1]: # same column, vertical conduction
        rwall = (domain.r[sink] + domain.r[source]) / 2
        ri = Q_(min(domain.r[source] if not sinkSide else rwall, domain.r[sink] if sinkSide else rwall), unitReg.inch).to(unitReg.foot)
        ro = Q_(max(domain.r[source] if not sinkSide else rwall, domain.r[sink] if sinkSide else rwall), unitReg.inch).to(unitReg.foot)
        return np.log(ro / ri) / (2 * np.pi * L * k)

def ConvectionHalfResistor(domain: domain_mmap.DomainMMAP, sink: tuple[int, int], source: tuple[int, int]) -> pint.Quantity:
    sinkSide = domain.material[sink] in MaterialType.COOLANT | MaterialType.EXHAUST
    convectionCell = sink if sinkSide else source
    isHoriz = source[0] == sink[0] # same row, horizontal convection
    sinkTop = sink[0] > source[0] # is the sink on top of the source
    mat = domain.material[convectionCell]

    if domain.material[convectionCell] in MaterialType.COOLANT:
        convectionCoeff = cooling2d.internal_flow_convection((domain.temperature[convectionCell]).to(unitReg.degR),(domain.pressure[convectionCell]).to(unitReg.psi), domain.area[convectionCell], domain.hydraulicDiameter[convectionCell])
        area = CoolantConvectionArea(domain, convectionCell[0], convectionCell[1], isHoriz, sinkTop, sinkSide)
    else:
        convectionCoeff = cooling2d.combustion_convection(domain.temperature[convectionCell].to(unitReg.degR), domain.velocity[convectionCell].to(unitReg.foot/unitReg.second))
        area = CombustionConvectionArea(domain, convectionCell[0], convectionCell[1], isHoriz, sinkTop, sinkSide)
    
    return 1 / (convectionCoeff * area)

def CoolantConvectionArea(domain: domain_mmap.DomainMMAP, row: int, col: int, isHoriz: bool, sinkTop: bool, sinkSide: bool) -> pint.Quantity: #TODO update
    wallPoint = (row, col)# if domain.material[row, col] in MaterialType.COOLANT_WALL else tuple(domain.previousFlow[row, col])
    innerRadius = Q_(domain.r[wallPoint] - domain.rstep/2, unitReg.inch).to(unitReg.foot)
    outerRadius = Q_(domain.r[wallPoint] + domain.rstep/2, unitReg.inch).to(unitReg.foot)

    coolingRatio = getCoolingCoverage(innerRadius, domain.coolingLoops[int(domain.id[row, col])])
    if isHoriz:
        return (outerRadius**2 - innerRadius**2) * coolingRatio
    else:
        outer = (sinkTop ^ sinkSide)               
        xstep = Q_(domain.xstep, unitReg.inch).to(unitReg.foot)           
        return (2 * np.pi * outerRadius * coolingRatio) * xstep if outer else (2 * np.pi * innerRadius * coolingRatio) * xstep

def CombustionConvectionArea(domain: domain_mmap.DomainMMAP, row: int, col: int, isHoriz: bool, sinkTop: bool, sinkSide: bool) -> pint.Quantity:
    innerRadius = Q_(domain.r[row, col] - domain.rstep/2, unitReg.inch).to(unitReg.foot)
    outerRadius = Q_(domain.r[row, col] + domain.rstep/2, unitReg.inch).to(unitReg.foot)
    if isHoriz:
        return np.pi * (outerRadius**2 - innerRadius**2)
    else:
        outer = (sinkTop ^ sinkSide)
        xstep = Q_(domain.xstep, unitReg.inch).to(unitReg.foot)
        return 2 * np.pi * outerRadius * xstep if outer else 2 * np.pi * innerRadius * xstep

def getCoolingCoverage(innerRadius: pint.Quantity, coolingLoopData: domain.CoolingLoop):
    numChannels = coolingLoopData.numChannels
    landWidth = coolingLoopData.landWidth
    return (2*np.pi*innerRadius/numChannels - landWidth)*numChannels/(2*np.pi*innerRadius)

def CombinationResistor(domain: domain_mmap.DomainMMAP, sink: tuple[int, int], source: tuple[int, int]) -> pint.Quantity:
    if domain.material[sink] == domain.material[source]:
        return ConductionResistor(domain, sink, source)

    if domain.material[source] in MaterialType.SOLID:
        sourceR = ConductionHalfResistor(domain, sink, source, False)
        # print("half solid")
    elif domain.material[source] in MaterialType.FLUID:
        sourceR = ConvectionHalfResistor(domain, sink, source)
        # print("half fluid")
    elif domain.material[source] in MaterialType.ADIABATIC:
        # print("half free")
        sourceR = Q_(2e31, unitReg.hour * unitReg.degR / unitReg.BTU)
    else:
        raise ValueError("Material not recognized")
    
    if domain.material[sink] in MaterialType.SOLID:
        sinkR = ConductionHalfResistor(domain, sink, source, True)
        # print("half solid")
    elif domain.material[sink] in MaterialType.FLUID:
        sinkR = ConvectionHalfResistor(domain, sink, source)
        # print("half fluid")
    elif domain.material[sink] in MaterialType.ADIABATIC:
        # print("half free")
        sinkR = Q_(2e31, unitReg.hour * unitReg.degR / unitReg.BTU)
    else:
        raise ValueError("Material not recognized")
    
    return (sourceR + sinkR).to(unitReg.hour * unitReg.degR / unitReg.BTU)

def CalculateCoolantPrimaryWall(domain: domain_mmap.DomainMMAP, row: int, col: int, solverSettings: dict, resistorSet: list[ThermalResistor] = []) -> list[CellUpdates]:
    wallCellUpdate = CellUpdates(row, col)
    prevFlow = tuple(domain.previousFlow[row, col])
    futureFlow = tuple(domain.futureFlow[row, col])
    coolingLoopData = domain.coolingLoops[int(domain.id[row, col])]
    if solverSettings.get("temperature", True) and domain.material[row, col] not in MaterialType.STATIC_TEMP:
        cp = fluid.get_cp(coolingLoopData.fluid, domain.temperature[row, col])
        mdotChannels = coolingLoopData.mdot
        numChannels = coolingLoopData.numChannels

        offsets = [(0, -1), (-1, 0), (1, 0), (0, 1)]
        for offset in offsets:
            if 0 <= row + offset[0] < domain.vpoints and 0 <= col + offset[1] < domain.hpoints:
                if domain.material[row + offset[0], col + offset[1]] in MaterialType.WALL:
                    resistorSet.append(ThermalResistor(CombinationResistor(domain, (row, col), (row + offset[0], col + offset[1])), domain.temperature[row + offset[0], col + offset[1]]))
        
        L1 = Q_(np.sqrt((domain.x[row, col] - domain.x[prevFlow])**2 + (domain.r[row, col] - domain.r[prevFlow])**2), unitReg.inch).to(unitReg.foot)
        L2 = Q_(np.sqrt((domain.x[row, col] - domain.x[futureFlow])**2 + (domain.r[row, col] - domain.r[futureFlow])**2), unitReg.inch).to(unitReg.foot)
        k = getConductivity(domain, row, col)

        resistorSet.append(ThermalResistor(L1/(domain.area[prevFlow].to(unitReg.foot**2) * numChannels * k), domain.temperature[prevFlow]))
        resistorSet.append(ThermalResistor(L2/(domain.area[futureFlow].to(unitReg.foot**2) * numChannels * k), domain.temperature[futureFlow]))

        num = sum([res.T / res.R for res in resistorSet])
        den = sum([1/res.R for res in resistorSet])

        # Tnew = ((num + mdotChannels*cp*domain.temperature[prevFlow])/(den + mdotChannels*cp)).to(unitReg.degR)
        # Tnew = ((num - .5*mdotChannels*cp*(domain.temperature[futureFlow] - domain.temperature[prevFlow])) / (den)).to(unitReg.degR)

        num = sum([(res.T - domain.temperature[row,col]) / res.R for res in resistorSet])
        Tnew = (num / (mdotChannels*cp)) + domain.temperature[prevFlow]
        wallCellUpdate.temperature = Tnew.to(domain.units["temperature"])

        if wallCellUpdate.temperature.m > 1e4 or wallCellUpdate.temperature.m < 0:
            print(domain.temperature[prevFlow])
            print(domain.temperature[futureFlow])
            print(mdotChannels*cp*(domain.temperature[futureFlow] - domain.temperature[prevFlow]))
            print(f"Temp out of bounds: {wallCellUpdate.temperature}")
            print(f"Row: {row}, Col: {col}")
            print(f"material: {domain.material[row, col]}")
            print(f"border: {domain.border[row, col]}")


    if solverSettings.get("pressure", True) and domain.material[row, col] not in MaterialType.STATIC_PRESS:
        futurePress = domain.pressure[futureFlow]
        deltaL = Q_(np.sqrt((domain.x[row, col] - domain.x[futureFlow])**2 + (domain.r[row, col] - domain.r[futureFlow])**2), unitReg.inch).to(unitReg.foot)

        (mu, cp, _, _, rho, _, _, _, _) = fluid.get_fluid_properties(coolingLoopData.fluid, domain.temperature[row, col], domain.pressure[row, col])
        mdotperchannel = coolingLoopData.mdot / coolingLoopData.numChannels
        Re = mdotperchannel*domain.hydraulicDiameter[row, col]/mu/domain.area[row, col]    # Reynolds number
        Re = Re.to(unitReg.dimensionless)
        f = fluid.DarcyFrictionFactor(Re, Q_(.05, unitReg.inch), domain.hydraulicDiameter[row, col])
        vel = (mdotperchannel / rho / domain.area[row, col]).to(unitReg.foot/unitReg.second)
        dp = fluid.FrictionPressureLoss(f, deltaL, domain.hydraulicDiameter[row, col], rho, vel)
        wallCellUpdate.pressure = (futurePress + dp).to(unitReg.psi)


    return [wallCellUpdate]

def CalculateCoolantBulkWall(domain: domain_mmap.DomainMMAP, row: int, col: int, solverSettings: dict) -> list[CellUpdates]:
    pass
    #steps:
    #    we enter the function with row, col pointing to a bulk flow that is on the wall and points to the series of other bulk on the same wall in future flow
    #    previous flow points to the associated wall point, which points to future and previous flows wall point. This wall point also stores channel area and id.
    # 1. traverse bulk walls and construct the resistor set for each convection point, conduction to other coolant is ignored (handled later)
    # 2. add wall flow to resistor set
    # 3. add conduction from coolant to resistor set, this uses flow area of previous and future flow in the wall point
    # 4. use coolant id to do mass flow transfers. We can assume a change in Cp is negligible, so we can use the same Cp for both in and out
    # 5. combine it all together
    # 
    # 2, 3, 4, 5 are needed to do wall coolant anyways, so we can just pass wall coolant a preloaded resistor set and it will do the rest
    #

    offsets = [(0, -1), (-1, 0), (1, 0), (0, 1)]

    resistorSet = []
    currentPoint = (row, col)
    while tuple(domain.futureFlow[currentPoint]) != currentPoint:
        for offset in offsets:
            if 0 <= currentPoint[0] + offset[0] < domain.vpoints and 0 <= currentPoint[1] + offset[1] < domain.hpoints:
                if domain.material[currentPoint[0] + offset[0], currentPoint[1] + offset[1]] in MaterialType.WALL:
                    resistorSet.append(ThermalResistor(CombinationResistor(domain, currentPoint, (currentPoint[0] + offset[0], currentPoint[1] + offset[1])), domain.temperature[currentPoint[0] + offset[0], currentPoint[1] + offset[1]]))
        currentPoint = tuple(domain.futureFlow[currentPoint])
    
    wallPoint = tuple(domain.previousFlow[row, col])
    cellUpdateArray = CalculateCoolantPrimaryWall(domain, wallPoint[0], wallPoint[1], solverSettings, resistorSet.copy())
    prevFlow = tuple(domain.previousFlow[row, col])
    cellUpdateArray.append(CellUpdates(row, col, temperature=domain.temperature[prevFlow], pressure=domain.pressure[prevFlow]))

    return cellUpdateArray

def CoolantBoundaryConditions(domain: domain_mmap.DomainMMAP, row: int, col: int, solverSettings: dict) -> list[CellUpdates]:
    offsets = [(0, -1), (-1, 0), (1, 0), (0, 1)]

    maxPress = 0
    maxTemp = 0

    for offset in offsets:
        if 0 <= row + offset[0] < domain.vpoints and 0 <= col + offset[1] < domain.hpoints:
            if domain.material[row + offset[0], col + offset[1]] in MaterialType.COOLANT:
                maxPress = max(maxPress, domain.pressure[row + offset[0], col + offset[1]].m)
                maxTemp = max(maxTemp, domain.temperature[row + offset[0], col + offset[1]].m)

    if domain.material[row,col] == DomainMaterial.COOLANT_INLET:
        return [CellUpdates(row, col, pressure=Q_(maxPress, unitReg.psi))]
    if domain.material[row,col] == DomainMaterial.COOLANT_OUTLET:
        return [CellUpdates(row, col, temperature=Q_(maxTemp, unitReg.degR))]

def CalculateBorderResistors(domain: domain_mmap.DomainMMAP, row: int, col: int, solverSettings: dict) -> list[ThermalResistor]:
    resSet = []

    offsets = [(0, -1), (-1, 0), (1, 0), (0, 1)]

    for offset in offsets:
        if 0 <= row + offset[0] < domain.vpoints and 0 <= col + offset[1] < domain.hpoints:
            R = CombinationResistor(domain, (row, col), (row + offset[0], col + offset[1]))
            resSet.append(ThermalResistor(R, domain.temperature[row + offset[0], col + offset[1]]))

    return resSet

def ConductionCoreResistors(domain: domain_mmap.DomainMMAP, row: int, col: int, solverSettings: dict) -> list[ThermalResistor]:
    resSet = []

    offsets = [(0, -1), (-1, 0), (1, 0), (0, 1)]

    for offset in offsets:
        if 0 <= row + offset[0] < domain.vpoints and 0 <= col + offset[1] < domain.hpoints:
            R = ConductionResistor(domain, (row, col), (row + offset[0], col + offset[1]))
            resSet.append(ThermalResistor(R, domain.temperature[row + offset[0], col + offset[1]]))

    return resSet

def MergeResistors(resSet: list[ThermalResistor]) -> pint.Quantity:
    if len(resSet) == 0:
        return Q_(-1, unitReg.degR)
    TRsum = Q_(0, str((resSet[0].T/resSet[0].R).units))
    Rsum = Q_(0, str(1/resSet[0].R.units))
    for res in resSet:
        if res.R.m > 1e31:
            continue
        TRsum += res.T / res.R
        Rsum += 1/res.R
    return (TRsum / Rsum).to(unitReg.degR)

def CalculateCell(domain: domain_mmap.DomainMMAP, row: int, col: int, **solverSettings: dict) -> list[CellUpdates]:
    if domain.material[row,col] in MaterialType.STATIC:
        return []
    
    if domain.r[row,col] - domain.rstep/2 <= 0:
        return []
    
    if domain.material[row,col] in MaterialType.WALL:
        if domain.border[row,col]:
            res = CalculateBorderResistors(domain, row, col, solverSettings)
            Tnew = MergeResistors(res)
            # if Tnew.m > 1e3:
            #     print(f"Wall temp out of bounds: {Tnew}")
            #     print(f"Row: {row}, Col: {col}")
            #     print(f"material: {domain.material[row, col]}")
            #     print(f"border: {domain.border[row, col]}")
            return [CellUpdates(row, col, temperature=Tnew)]
        res = ConductionCoreResistors(domain, row, col, solverSettings)
        Tnew = MergeResistors(res)
        # if Tnew.m > 1e3:
        #     print(f"Wall temp out of bounds: {Tnew}")
        #     print(f"Row: {row}, Col: {col}")
        #     print(f"material: {domain.material[row, col]}")
        #     print(f"border: {domain.border[row, col]}")
        return [CellUpdates(row, col, temperature=Tnew)]
    
    if domain.material[row,col] in MaterialType.COOLANT:
        if domain.material[row,col] == DomainMaterial.COOLANT_BULK:
            if domain.border[row,col]:
                return CalculateCoolantBulkWall(domain, row, col, solverSettings) #TODO implement
            prevFlow = tuple(domain.previousFlow[row, col])
            return [CellUpdates(row, col, temperature=domain.temperature[prevFlow], pressure=domain.pressure[prevFlow])]
        if domain.material[row,col] == DomainMaterial.COOLANT_WALL:
            return CalculateCoolantPrimaryWall(domain, row, col, solverSettings) #TODO implement
        if domain.material[row,col] in {DomainMaterial.COOLANT_INLET, DomainMaterial.COOLANT_OUTLET}:
            return CoolantBoundaryConditions(domain, row, col, solverSettings)
    raise ValueError(f"Material not recognized {domain.material[row, col]}")