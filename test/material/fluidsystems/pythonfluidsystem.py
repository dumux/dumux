from dumux.material import Component

# Air = Component("Air")
# print(Air.name())
# print("molar mass", Air.molarMass())
# print("gas thermal conductivity", Air.gasThermalConductivity(300.0, 1e5))
#
#
# Calcite = Component("Calcite")
# print(Calcite.name())
# print("molar mass", Calcite.molarMass())
# print("solid density", Calcite.solidDensity(300.0))
#
#
H2O = Component("H2O")
# print(H2O.name())
# print("molar mass", H2O.molarMass())
# print("vapor pressure", H2O.vaporPressure(300.0))


def density(temperature, pressure):
    return 1e-5*pressure

def diffusionCoefficient(fluidstate, phaseIdx, compIdx):
    return 1e-9

def enthalpy(temperature, pressure):
    return H2O.liquidEnthalpy(temperature, pressure)
