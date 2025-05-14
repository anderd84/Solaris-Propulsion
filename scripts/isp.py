from ambiance import Atmosphere
import matplotlib.pyplot as plt
from icecream import ic
import numpy as np
import joblib

from general.units import Q_, unitReg
from general.units import PSI, FT, IN, S, LBF, LBM
import general.design as DESIGN
from nozzle import plug
from nozzle import analysis

def getIsp(i, pamb, cont, exhaust, Tt, Rt, Re, mdot):
    from general.units import Q_, unitReg
    from general.units import PSI, FT, IN, S, LBF, LBM
    g0 = Q_(32.2, FT/S/S)
    janusVacIsp = Q_(304, S)
    janusmdot = Q_(7.5, LBM/S)
    janusVacThrust = janusVacIsp * janusmdot * g0
    janusAe = np.pi * (Q_(6.25, IN)/2)**2
    exhaust = DESIGN.exhaustGas

    Tt = Q_(Tt, unitReg.radian)
    Rt = Q_(Rt, unitReg.inch)
    Re = Q_(Re, unitReg.inch)
    mdot = Q_(mdot, LBM/S)

    pamb = Q_(pamb, PSI)
    rlines, llines, streams = analysis.CalculateComplexField(cont, pamb, exhaust, 1, Tt, Rt, Re.magnitude, 50, 0, 3 if pamb.magnitude > 6.75 else 1)
    T = analysis.CalculateThrust(exhaust, pamb, Tt, Rt, Re, streams[0], cont[-1].r)
    print(f"Spike {i}: {T.to(LBF)}")
    print(f"Janus {i}: {(janusVacThrust - pamb*janusAe).to(LBF)}")

    ispSpike = ((T/mdot)/g0).to(S).magnitude

    ispJanus = (((janusVacThrust - pamb*janusAe)/janusmdot)/g0).to(S).magnitude
    return (i, ispSpike, ispJanus)

def main():
    Re = Q_(3.2, IN)
    janusVacThrust = Q_(2178, unitReg.pound_force)
    janusmdot = Q_(7.5, unitReg.pound/unitReg.seconds)
    janusAe = Q_(30.87, unitReg.inch**2)
    exhaust = DESIGN.exhaustGas

    cont, field, outputData = plug.CreateRaoContour(exhaust, DESIGN.chamberPressure, DESIGN.designAmbientPressure, DESIGN.basePressure, Re, DESIGN.lengthMax)
    Rt = outputData["radiusThroat"]
    Tt = outputData["thetaThroat"]
    Re = outputData["radiusLip"]
    ic(outputData["areaRatio"])

    phi = np.pi/2 + Tt
    Astar = np.pi/np.sin(phi) * (Re**2 - Rt**2)

    h0 = Q_(3500, FT)
    hf = Q_(30000, FT)

    mdot = DESIGN.totalmdot

    pts = 30
    ispSpike = np.zeros(pts)
    ispJanus = np.zeros(pts)

    alts = [min(h.to(unitReg.meter).magnitude, 81000) for h in np.linspace(h0, hf, pts)]
    pressures = [Q_(p, unitReg.pascal).to(unitReg.psi).magnitude for p in Atmosphere(alts).pressure]

    with joblib.Parallel(n_jobs=15) as parallel:

        outputs = parallel(joblib.delayed(getIsp)(i, p, cont, exhaust, Tt.magnitude, Rt.magnitude, Re.magnitude, mdot.magnitude) for i, p in enumerate(pressures))

        for i, _ispSpike, _ispJanus in outputs:
            ic(_ispSpike, _ispJanus)
            ispSpike[i] = _ispSpike
            ispJanus[i] = _ispJanus

    # ic(isp)
    plt.plot(np.array(alts)*3.28084, ispSpike, '-r', linewidth=2)
    plt.plot(np.array(alts)*3.28084, ispJanus, '-b', linewidth=2)
    plt.xlabel("Altitude (ft)")
    plt.ylabel("ISP (s)")
    plt.legend(["Aerospike", "Janus 4.2"])
    ic(np.mean(ispSpike), np.mean(ispJanus))
    plt.show()


if __name__ == "__main__":
    main()