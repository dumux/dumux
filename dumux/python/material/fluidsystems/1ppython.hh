#ifndef DUMUX_PYTHON_MATERIAL_FLUIDSYSTEM_HH
#define DUMUX_PYTHON_MATERIAL_FLUIDSYSTEM_HH

#include <string>
#include <memory>
#include <tuple>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/componenttraits.hh>
#include <dumux/python/material/fluidstate.hh>
#include <dumux/io/name.hh>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/embed.h>

#include <dumux/python/material/fluidstate.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux::Python::FluidSystems {

namespace Detail {
template <std::size_t size, class F, std::size_t I = 0>
inline auto switchSequenceCall([[maybe_unused]] std::size_t i, F&& f)
-> decltype(f(std::integral_constant<std::size_t, 0>{}))
{
    if constexpr(I < size) {
        switch (i) {
            case I: return f(std::integral_constant<std::size_t, I>{});
            default: return switchSequenceCall<size, F, I+1>(i, std::forward<F>(f));
        }
    }
    else
        DUNE_THROW(Dune::InvalidStateException, "Index out of range!");
}

}
/*!
 * \ingroup Fluidsystems
 * \brief A liquid phase consisting of a single component
 */
template <class Scalar, class PythonFluidSystemName,  class FirstComponent, class ...OtherComponents>
class OnePLiquid
: public Dumux::FluidSystems::Base<Scalar, OnePLiquid<Scalar, PythonFluidSystemName, FirstComponent, OtherComponents...> >
{
    using ThisType = OnePLiquid<Scalar, PythonFluidSystemName, FirstComponent, OtherComponents...>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

    static_assert(ComponentTraits<FirstComponent>::hasLiquidState);

    // static_assert(std::conjunction_v<std::integral_constant<bool, ComponentTraits<OtherComponents>::hasLiquidState>...>,
    //               "One of the components does not implement a liquid state!");

public:
    template<std::size_t i>
    using Component = typename std::tuple_element_t<i, std::tuple<FirstComponent, OtherComponents...>>;

    using ParameterCache = NullParameterCache;

    static constexpr std::size_t numPhases = 1;  //!< Number of phases in the fluid system
    static constexpr std::size_t numComponents = sizeof...(OtherComponents) + 1; //!< Number of components in the fluid system

    static constexpr std::size_t phase0Idx = 0; //!< index of the only phase
    static constexpr std::size_t comp0Idx = 0; //!< index of the only component

    /*!
     * \brief Initialize the fluid system's static parameters generically
     */
    static void init()
    { }

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return IOName::liquidPhase(); }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx = 0)
    {
        return Detail::switchSequenceCall<numComponents>(compIdx, [](auto _compIdx)
        {
            constexpr auto compIdx = std::decay_t<decltype(_compIdx)>::value;
            return Component<compIdx>::name();
        });
    }

    /*!
     * \brief There is only one phase, so not mass transfer between phases can occur
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief Returns whether the fluid is a liquid
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return false; }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if only a single component is involved. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealMixture(int phaseIdx = 0)
    { return true; }

    /*!
     * \brief Returns true if the fluid is assumed to be compressible
     */
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return FirstComponent::liquidIsCompressible(); }
    // { return std::disjunction_v<std::integral_constant<bool, ComponentT::liquidIsCompressible()>...>; }

    /*!
     * \brief Returns true if the fluid viscosity is constant
     */
    static constexpr bool viscosityIsConstant(int phaseIdx = 0)
    { return FirstComponent::liquidViscosityIsConstant(); }
    // { return std::conjunction_v<std::integral_constant<bool, ComponentT::liquidViscosityIsConstant()>...>; }

    /*!
     * \brief Returns true if the fluid is assumed to be an ideal gas
     */
    static constexpr bool isIdealGas(int phaseIdx = 0)
    { return false; /* we're a liquid! */ }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass(int compIdx = 0)
    {
        return Detail::switchSequenceCall<numComponents>(compIdx, [](auto _compIdx)
        {
            constexpr auto compIdx = std::decay_t<decltype(_compIdx)>::value;
            return Component<compIdx>::molarMass();
        });
    }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of the component
     */
    static Scalar criticalTemperature(int compIdx = 0)
    {
        return Detail::switchSequenceCall<numComponents>(compIdx, [](auto _compIdx)
        {
            constexpr auto compIdx = std::decay_t<decltype(_compIdx)>::value;
            return Component<compIdx>::criticalTemperature();
        });
    }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of the component
     */
    static Scalar criticalPressure(int compIdx = 0)
    {
        return Detail::switchSequenceCall<numComponents>(compIdx, [](auto _compIdx)
        {
            constexpr auto compIdx = std::decay_t<decltype(_compIdx)>::value;
            return Component<compIdx>::criticalPressure();
        });
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature(int compIdx = 0)
    {
        return Detail::switchSequenceCall<numComponents>(compIdx, [](auto _compIdx)
        {
            constexpr auto compIdx = std::decay_t<decltype(_compIdx)>::value;
            return Component<compIdx>::tripleTemperature();
        });
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure(int compIdx = 0)
    {
        return Detail::switchSequenceCall<numComponents>(compIdx, [](auto _compIdx)
        {
            constexpr auto compIdx = std::decay_t<decltype(_compIdx)>::value;
            return Component<compIdx>::triplePressure();
        });
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        return pythonFluidSytem_().attr("vaporPressure")(temperature).template cast<Scalar>();
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("density")(temperature, pressure).template cast<Scalar>();
    }

    using Base::density;
    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    template <class FluidState>
    static Scalar density(const FluidState& fluidState,
                          const int phaseIdx)
    {
        return pythonFluidSytem_().attr("density")(fluidState, phaseIdx).template cast<Scalar>();
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the main component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, const int phaseIdx)
    {
        return pythonFluidSytem_().attr("molarDensity")(fluidState, phaseIdx).template cast<Scalar>();
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    static Scalar molarDensity(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("molarDensity")(temperature, pressure).template cast<Scalar>();
    }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of the component at a given density and temperature.
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {
        return pythonFluidSytem_().attr("pressure")(temperature, pressure).template cast<Scalar>();
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("enthalpy")(temperature, pressure).template cast<Scalar>();
    }

    using Base::enthalpy;
    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const int phaseIdx)
    {
        return enthalpy(fluidState.temperature(phaseIdx),
                        fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("internalEnergy")(temperature, pressure).template cast<Scalar>();
    }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("viscosity")(temperature, pressure).template cast<Scalar>();
    }

    using Base::viscosity;
    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const int phaseIdx)
    {
        return pythonFluidSytem_().attr("viscosity")(fluidState, phaseIdx).template cast<Scalar>();
    }

    using Base::fugacityCoefficient;
    /*!
     * \copybrief Base::fugacityCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    using Base::diffusionCoefficient;
    /*!
     * \copybrief Base::diffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        return pythonFluidSytem_().attr("diffusionCoefficient")(fluidState, phaseIdx, compIdx).template cast<Scalar>();
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \copybrief Base::binaryDiffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the component to consider
     * \param compJIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        return pythonFluidSytem_().attr("binaryDiffusionCoefficient")(fluidState, phaseIdx, compIIdx, compJIdx).template cast<Scalar>();
    }

    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("thermalConductivity")(temperature, pressure).template cast<Scalar>();
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        return pythonFluidSytem_().attr("thermalConductivity")(fluidState, phaseIdx).template cast<Scalar>();
    }

    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    {
        return pythonFluidSytem_().attr("heatCapacity")(temperature, pressure).template cast<Scalar>();
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const int phaseIdx)
    {
        return heatCapacity(fluidState.temperature(phaseIdx),
                            fluidState.pressure(phaseIdx));
    }
private:

    static pybind11::object& pythonFluidSytem_()
    {
        static pybind11::scoped_interpreter guard{};
        static pybind11::object pyFluidSystem = pybind11::module::import(PythonFluidSystemName::get());

        [[maybe_unused]] static const bool registered = [&]()
        {
            using FS = Dumux::CompositionalFluidState<Scalar, ThisType>;
            Dumux::Python::Impl::registerFluidState<FS>(pyFluidSystem);
            return true;
        }();

        return pyFluidSystem;
    }
};

}


#endif
