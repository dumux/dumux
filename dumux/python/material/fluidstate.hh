#ifndef DUMUX_PYTHON_MATERIAL_FLUIDSTATE_HH
#define DUMUX_PYTHON_MATERIAL_FLUIDSTATE_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dune/python/common/typeregistry.hh>
#include <dune/common/classname.hh>

namespace Dumux::Python {

namespace Impl {

template<class FluidState>
void registerFluidState(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<FluidState>(
        scope, "FluidState",
        GenerateTypeName(Dune::className<FluidState>()),
        IncludeFiles{}
    );

    if (addedToRegistry)
    {
        using pybind11::operator""_a;

        cls.def("temperature", [](FluidState& self, const int phaseIdx) {
            return self.temperature(0);
        });

        cls.def("temperature", [](FluidState& self) {
            return self.temperature();
        });

        cls.def("fugacity", [](FluidState& self, const int compIdx) {
            return self.fugacity(compIdx);
        });

        cls.def("fugacity", [](FluidState& self, const int phaseIdx, int compIdx) {
            return self.fugacity(phaseIdx, compIdx);
        });

        cls.def("wettingPhase", &FluidState::wettingPhase);
        cls.def("pressure", &FluidState::pressure, "phaseIdx"_a);
        cls.def("moleFraction", &FluidState::moleFraction, "phaseIdx"_a, "compIdx"_a);
        cls.def("massFraction", &FluidState::massFraction, "phaseIdx"_a, "compIdx"_a);
        cls.def("molarity", &FluidState::molarity, "phaseIdx"_a, "compIdx"_a);
        cls.def("fugacityCoefficient", &FluidState::fugacityCoefficient, "phaseIdx"_a, "compIdx"_a);
        cls.def("partialPressure", &FluidState::partialPressure, "phaseIdx"_a, "compIdx"_a);
        cls.def("molarVolume", &FluidState::molarVolume, "phaseIdx"_a);
        cls.def("density", &FluidState::density, "phaseIdx"_a);
        cls.def("molarDensity", &FluidState::molarDensity, "phaseIdx"_a);
        cls.def("averageMolarMass", &FluidState::averageMolarMass, "phaseIdx"_a);
        cls.def("phaseMassFraction", &FluidState::phaseMassFraction, "phaseIdx"_a);
        cls.def("saturation", &FluidState::saturation, "phaseIdx"_a);
        cls.def("enthalpy", &FluidState::enthalpy, "phaseIdx"_a);
        cls.def("internalEnergy", &FluidState::internalEnergy, "phaseIdx"_a);
        cls.def("viscosity", &FluidState::viscosity, "phaseIdx"_a);
    }
}
}

} // namespace Dumux::Python

#endif
