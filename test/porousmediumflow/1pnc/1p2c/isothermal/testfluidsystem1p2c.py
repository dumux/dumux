from dumux.material import Component

liquidPhaseIdx = 0
H2OIdx = 0
N2Idx = 1


H2O = Component("H2O")
N2 = Component("N2")



def density(fluidState, phaseIdx):
    T = fluidState.temperature(phaseIdx)
    p = fluidState.pressure(phaseIdx)
    return H2O.liquidMolarDensity(T, p) * (H2O.molarMass()*fluidState.moleFraction(liquidPhaseIdx, H2OIdx) + N2.molarMass()*fluidState.moleFraction(liquidPhaseIdx, N2Idx))

def molarDensity(fluidState, phaseIdx):
    T = fluidState.temperature(phaseIdx)
    p = fluidState.pressure(phaseIdx)
    return H2O.liquidMolarDensity(T, p);

def viscosity(fluidState, phaseIdx):
    T = fluidState.temperature(phaseIdx)
    p = fluidState.pressure(phaseIdx)
    return H2O.liquidViscosity(T, p)

def binaryDiffusionCoefficient(fluidState, phaseIdx, compIIdx, compJIdx):
    if compIIdx > compJIdx:
        swap(compIIdx, compJIdx);

    Texp = 273.15 + 25
    Dexp = 2.01e-9
    return Dexp * fluidState.temperature(phaseIdx)/Texp;
