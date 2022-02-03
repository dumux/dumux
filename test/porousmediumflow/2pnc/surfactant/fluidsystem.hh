#ifndef DUMUX_TEST_2P3C_SURFACTANT_FLUIDSYSTEM_HH
#define DUMUX_TEST_2P3C_SURFACTANT_FLUIDSYSTEM_HH

#include <cassert>

#include <dumux/material/fluidsystems/base.hh>

namespace Dumux::FluidSystems {

template <class Scalar>
class TestSurfactant
: public Base<Scalar, TestSurfactant<Scalar>>
{
    using ThisType = TestSurfactant<Scalar>;
    using Base = FluidSystems::Base<Scalar, ThisType>;

public:
    // Two phases: wetting (0) and oil/nonwetting (1)
    static constexpr int numPhases = 2;

    static constexpr int wettingPhaseIdx = 0;
    static constexpr int nonwettingPhaseIdx = 1;

    static constexpr int phase0Idx = 0;
    static constexpr int phase1Idx = 1;

    static constexpr int numComponents = 3;

    static constexpr int waterCompIdx = 0;
    static constexpr int oilCompIdx = 1;
    static constexpr int surfactantCompIdx = 2;

    static constexpr int comp0Idx = 0;
    static constexpr int comp1Idx = 1;
    static constexpr int comp2Idx = 2;


    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == wettingPhaseIdx)
            return "w";
        return "n";
    }

    static constexpr bool isMiscible()
    { return false; }

    static constexpr bool isGas(int phaseIdx)
    { return false; }

    static constexpr bool isIdealMixture(int phaseIdx)
    { return true; }

    static constexpr bool isCompressible(int phaseIdx)
    { return false; }

    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return true; }

    /****************************************
     * Component related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == waterCompIdx)
            return "water";
        if (compIdx == surfactantCompIdx)
            return "surfactant";

        return "oil";
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == waterCompIdx)
            return 18e-3;
        if (compIdx == surfactantCompIdx)
            return 18e-3;

        // oil component
        return 0.350;
    }

    using Base::density;
    /*!
     * \brief Calculate the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     */
    template <class FluidState>
    static Scalar density(const FluidState& fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == wettingPhaseIdx)
        {
            static const Scalar wettingPhaseDensity = getParam<Scalar>("FluidSystem.WettingPhaseDensity");
            return wettingPhaseDensity;
        }

        static const Scalar nonwettingPhaseDensity = getParam<Scalar>("FluidSystem.NonwettingPhaseDensity");
        return nonwettingPhaseDensity;
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == wettingPhaseIdx)
            return density(fluidState, phaseIdx) / molarMass(waterCompIdx);

        return density(fluidState, phaseIdx) / molarMass(oilCompIdx);
    }

    using Base::viscosity;
    /*!
     * \brief Return the viscosity of a phase \f$\mathrm{[Pa*s]}\f$.
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        if (phaseIdx == wettingPhaseIdx)
        {
            static const Scalar wettingPhaseViscosity
                = getParam<Scalar>("FluidSystem.WettingPhaseViscosity");
            return wettingPhaseViscosity;
        }

        static const Scalar nonwettingPhaseViscosity
            = getParam<Scalar>("FluidSystem.NonwettingPhaseViscosity");
        return nonwettingPhaseViscosity;
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\mathrm{\phi^\kappa_\alpha}\f$ is connected to the
     * fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction \f$\mathrm{x^\kappa_\alpha}\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     * \f]
     *
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     * \param compIdx index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == wettingPhaseIdx)
            if (compIdx == waterCompIdx || compIdx == surfactantCompIdx)
                return 0.0; // water and solubles stay in water phase

        if (phaseIdx == nonwettingPhaseIdx)
            if (compIdx == oilCompIdx)
                return 0.0; // oil stays in oil phase

        return 100.0; // arbitrary?
    }

    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState, int phaseIdx, int compIIdx, int compJIdx)
    { return 0.0; }
};

} // end namespace Dumux::FluidSystems

#endif
