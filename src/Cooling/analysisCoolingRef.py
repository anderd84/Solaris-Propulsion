from enum import IntEnum
from Cooling import cooling2d, domain
import General.design as DESIGN
# from General.units import Q_, unitReg
from Cooling import cooling2d as cooling_func
from Cooling.material import DomainMaterial, MaterialType
import numpy as np
from icecream import ic
from General.units import Q_, Direction
import pint

unitReg = pint.get_application_registry()
Q_ = unitReg.Quantity

class ResistorSet:
    R: list[Q_] = [Q_(0, unitReg.hour * unitReg.degR / unitReg.BTU) for _ in range(4)]
    T: list[Q_] = [Q_(0, unitReg.degR) for _ in range(4)]

    def getSums(self):
        TRsum = Q_(0, str((self.T[0]/self.R[0]).units))
        Rsum = Q_(0, str(1/self.R[0].units))
        for i in range(4):
            if self.R[i] > Q_(2e8, unitReg.hour * unitReg.degR / unitReg.BTU) or self.R[i] <= Q_(2e-31, unitReg.hour * unitReg.degR / unitReg.BTU):
                continue
            TRsum += self.T[i] / self.R[i]
            Rsum += 1/self.R[i]
            # print(TRsum, Rsum)
        return TRsum, Rsum

    def getTnew(self):
        TRsum, Rsum = self.getSums()
        return TRsum / Rsum

    def getQin(self, Tprev: Q_):
        Qsum = Q_(0, str((self.T[0]/self.R[0]).units))
        for i in range(4):
            if self.R[i] > Q_(2e8, unitReg.hour * unitReg.degR / unitReg.BTU) or self.R[i] <= Q_(2e-31, unitReg.hour * unitReg.degR / unitReg.BTU):
                continue
            Qsum += (self.T[i] - Tprev)/self.R[i]
        return Qsum


def getConductivity(domain: domain.DomainMMAP, row: int, col: int):
    if domain.material[row,col] in MaterialType.WALL:
        return cooling_func.conduction_grcop(domain.temperature[row,col].to(unitReg.degR))
    else:
        return cooling_func.conduction_rp1(domain.temperature[row,col].to(unitReg.degR))

def ConductionResistor(domain: domain.DomainMMAP, sink: tuple[int, int], source: tuple[int, int]): # sink is to, source is from
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

def ConductionHalfResistor(domain: domain.DomainMMAP, sink: tuple[int, int], source: tuple[int, int], sinkSide: bool = True):
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

def ConvectionHalfResistor(domain: domain.DomainMMAP, sink: tuple[int, int], source: tuple[int, int]):
    sinkSide = domain.material[sink] in MaterialType.COOLANT | MaterialType.EXHAUST
    convectionCell = sink if sinkSide else source
    isHoriz = source[0] == sink[0] # same row, horizontal convection
    sinkTop = sink[0] > source[0] # is the sink on top of the source
    mat = domain.material[convectionCell]

    if domain.material[convectionCell] in MaterialType.COOLANT:
        convectionCoeff = cooling_func.internal_flow_convection((domain.temperature[convectionCell]).to(unitReg.degR),(domain.pressure[convectionCell]).to(unitReg.psi), domain.area[convectionCell], domain.hydraulicDiameter[convectionCell])
        area = CoolantConvectionArea(domain, convectionCell[0], convectionCell[1], isHoriz, sinkTop, sinkSide)
        # print(convectionCoeff, area)
        # print(domain.flowHeight[convectionCell])
    else:
        convectionCoeff = cooling_func.combustion_convection(domain.temperature[convectionCell].to(unitReg.degR), domain.velocity[convectionCell].to(unitReg.foot/unitReg.second))
        area = CombustionConvectionArea(domain, convectionCell[0], convectionCell[1], isHoriz, sinkTop, sinkSide)
    
    return 1 / (convectionCoeff * area)

def CoolantConvectionArea(domain: domain.DomainMMAP, row: int, col: int, isHoriz: bool, sinkTop: bool, sinkSide: bool):
    wallPoint = (row, col)# if domain.material[row, col] in MaterialType.COOLANT_WALL else tuple(domain.previousFlow[row, col])
    innerRadius = Q_(domain.r[wallPoint] - domain.rstep/2, unitReg.inch).to(unitReg.foot)
    outerRadius = Q_(domain.r[wallPoint] + domain.rstep/2, unitReg.inch).to(unitReg.foot)
    landRadius = innerRadius

    theta = getCoolingArcAngle(landRadius)
    if isHoriz:
        return (outerRadius**2 - innerRadius**2) * theta/2 * DESIGN.NumberofChannels
    else:
        outer = (sinkTop ^ sinkSide)               
        xstep = Q_(domain.xstep, unitReg.inch).to(unitReg.foot)           
        return (2 * np.pi * outerRadius * theta * DESIGN.NumberofChannels) * xstep if outer else (2 * np.pi * innerRadius * theta * DESIGN.NumberofChannels) * xstep

def CombustionConvectionArea(domain: domain.DomainMMAP, row: int, col: int, isHoriz: bool, sinkTop: bool, sinkSide: bool):
    innerRadius = Q_(domain.r[row, col] - domain.rstep/2, unitReg.inch).to(unitReg.foot)
    outerRadius = Q_(domain.r[row, col] + domain.rstep/2, unitReg.inch).to(unitReg.foot)
    if isHoriz:
        return np.pi * (outerRadius**2 - innerRadius**2)
    else:
        outer = (sinkTop ^ sinkSide)
        xstep = Q_(domain.xstep, unitReg.inch).to(unitReg.foot)
        return 2 * np.pi * outerRadius * xstep if outer else 2 * np.pi * innerRadius * xstep
    
def getCoolingArcAngle(innerRadius: Q_):
    return 2*np.pi/DESIGN.NumberofChannels - (DESIGN.landWidth / innerRadius)

def CalculateCoreResistors(domain: domain.DomainMMAP, row: int, col: int):
    resSet = ResistorSet()

    if col > 0: # left
        resSet.R[Direction.LEFT] = ConductionResistor(domain, (row, col), (row, col-1))
        resSet.T[Direction.LEFT] = domain.temperature[row, col-1]

    if row > 0: # upper
        resSet.R[Direction.UPPER] = ConductionResistor(domain, (row, col), (row-1, col))
        resSet.T[Direction.UPPER] = domain.temperature[row-1, col]

    if row < domain.vpoints - 1: # bottom
        resSet.R[Direction.LOWER] = ConductionResistor(domain, (row, col), (row+1, col))
        resSet.T[Direction.LOWER] = domain.temperature[row+1, col]

    if col < domain.hpoints - 1: # right
        resSet.R[Direction.RIGHT] = ConductionResistor(domain, (row, col), (row, col+1))
        resSet.T[Direction.RIGHT] = domain.temperature[row, col+1]

    return resSet

def CombineResistors(resSet: ResistorSet):
    for i in range(4):
        if resSet.R[i] > Q_(1e10, unitReg.hour * unitReg.degR / unitReg.BTU):
            resSet.T[i] = Q_(0, unitReg.degR)
            resSet.R[i] = Q_(1, unitReg.hour * unitReg.degR / unitReg.BTU)
    return resSet.getTnew()

def GetResistor(domain: domain.DomainMMAP, sink: tuple[int, int], source: tuple[int, int]):
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
    
    return sourceR + sinkR

def CalculateBorderResistors(domain: domain.DomainMMAP, row: int, col: int):
    resSet = ResistorSet()

    if col > 0: # left
        resSet.R[Direction.LEFT] = GetResistor(domain, (row, col), (row, col-1))
        resSet.T[Direction.LEFT] = domain.temperature[row, col-1]

    if row > 0: # upper
        resSet.R[Direction.UPPER] = GetResistor(domain, (row, col), (row-1, col))
        resSet.T[Direction.UPPER] = domain.temperature[row-1, col]

    if row < domain.vpoints - 1: # bottom
        resSet.R[Direction.LOWER] = GetResistor(domain, (row, col), (row+1, col))
        resSet.T[Direction.LOWER] = domain.temperature[row+1, col]

    if col < domain.hpoints - 1: # right
        resSet.R[Direction.RIGHT] = GetResistor(domain, (row, col), (row, col+1))
        resSet.T[Direction.RIGHT] = domain.temperature[row, col+1]

    # print(resSet.R, resSet.T)
    return resSet

def CalculateCoolant(domain: domain.DomainMMAP, row: int, col: int):
    previousFlow = tuple(domain.previousFlow[row, col])

    resSet = CalculateBorderResistors(domain, row, col)
    Tprev = domain.temperature[row, col]

    deltaL = Q_(np.sqrt((domain.x[row, col] - domain.x[previousFlow])**2 + (domain.r[row, col] - domain.r[previousFlow])**2), unitReg.inch).to(unitReg.foot)
    flowHeight = (domain.flowHeight[row, col] + domain.flowHeight[previousFlow]) / 2

    return cooling2d.heatcoolant(domain.temperature[previousFlow], Tprev, resSet, domain.pressure[previousFlow], domain.pressure[row, col], domain.area[row, col], domain.hydraulicDiameter[row, col], deltaL)

def CalculateCell(domain: domain.DomainMMAP, row: int, col: int):
    if domain.material[row,col] in MaterialType.STATIC_TEMP:
        return domain.temperature[row,col]
    
    if domain.r[row,col] - domain.rstep/2 <= 0:
        return domain.temperature[row,col]
    
    if domain.material[row,col] in MaterialType.WALL:
        if not domain.border[row,col]:
            return CalculateCoreResistors(domain, row, col).getTnew()
        return CalculateBorderResistors(domain, row, col).getTnew()
    
    if domain.material[row,col] in MaterialType.COOLANT:
        if domain.material[row,col] == DomainMaterial.COOLANT_BULK:
            prevFlow = tuple(domain.previousFlow[row, col])
            return (domain.temperature[prevFlow], domain.pressure[prevFlow])
        return CalculateCoolant(domain, row, col)
    raise ValueError("Material not recognized")

def Cell(d: domain.DomainMMAP, row: int, col: int):
    out = CalculateCell(d, row, col)
    if isinstance(out, tuple):
        return out
    return (out, d.pressure[row, col])