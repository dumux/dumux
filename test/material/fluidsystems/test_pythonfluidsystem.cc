#include <config.h>

#include <dumux/python/material/fluidsystems/1ppython.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/embed.h>
#include <dumux/python/material/fluidstate.hh>

int main()
{
    struct Name
    {
        static constexpr auto get()
        { return "pythonfluidsystem"; }
    };
    using FS = Dumux::Python::FluidSystems::OnePLiquid<double, Name, Dumux::Components::Constant<0, double>, Dumux::Components::H2O<double>>;
    std::cout << "density is " << FS::density(300.0, 1e5) << std::endl;
    std::cout << "comp name 0 is " << FS::componentName(0) << std::endl;
    std::cout << "comp name 1 is " << FS::componentName(1) << std::endl;
    std::cout << "is compressible " << FS::isCompressible() << std::endl;

    Dumux::CompositionalFluidState<double, FS> fluidState;
    fluidState.setTemperature(300.0);
    fluidState.setPressure(0, 1e5);

    std::cout << "diffusion coefficient " << FS::diffusionCoefficient(fluidState, 0, 0) << std::endl;
    std::cout << "enthalpy " << FS::enthalpy(fluidState, 0) << std::endl;

    return 0;
}
