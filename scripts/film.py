from fluids import fluid
from general.units import Q_, unitReg
import general.design as DESIGN


hydroD = Q_(0.034, unitReg.inch)
area = Q_(.00125, unitReg.inch**2)  # area of one channel
(mu, cp, _, _, rho, _, _, _, _) = fluid.get_fluid_properties('RP-1', Q_(700, unitReg.degR), Q_(500, unitReg.psi))
mdotperchannel = DESIGN.Fuel_Total / 300
Re = mdotperchannel*hydroD/mu/area    # Reynolds number
Re = Re.to(unitReg.dimensionless)
f = fluid.DarcyFrictionFactor(Re, Q_(.05, unitReg.inch), hydroD)
vel = (mdotperchannel / rho / area).to(unitReg.foot/unitReg.second)
dp = fluid.FrictionPressureLoss(f, Q_(8, unitReg.inch), hydroD, rho, vel)
print(dp.to(unitReg.psi))
print(f)