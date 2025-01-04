import pint
import joblib

def main():
    a = pint.get_application_registry()
    print(dir(a))
    unitReg = pint.UnitRegistry()
    unitReg.default_system = 'US'
    unitReg.formatter.default_format = "~P"  # Compact unit formatting
    unitReg.define("lbmol = 453.59237 * mol")  # 1 lbmol = 453.59237 mol (since 1 lb = 453.59237 g)
    pint.set_application_registry(unitReg)

    print(hasattr(a, 'lbmol')) 


    with joblib.Parallel(3) as p:
        outs = p(joblib.delayed(testP)() for _ in range(5))

    sums = a.Quantity(0, 'inch')
    for out in outs:
        print(out)
        sums += out

    print(sums)


def testP():
    a = pint.get_application_registry()
    print(f"lb: {hasattr(a, 'lb')}")
    print(f"lbmol: {hasattr(a, 'lbmol')}")
    if not hasattr(a, 'lbmol'):
        print("Creating new registry")
        unitReg = pint.UnitRegistry()
        unitReg.default_system = 'US'
        unitReg.formatter.default_format = "~P"  # Compact unit formatting
        unitReg.define("lbmol = 453.59237 * mol")  # 1 lbmol = 453.59237 mol (since 1 lb = 453.59237 g)
        pint.set_application_registry(unitReg)
    print(a)
    return a.Quantity(1, 'inch')

if __name__ == "__main__":
    main()
