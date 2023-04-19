// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_COMPOSITIONAL_FLUID_STATE_HH
#define DUMUX_COMPOSITIONAL_FLUID_STATE_HH

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <cassert>
#include <array>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
template <class ScalarType, class FluidSystem>
class CompositionalFluidState
{
public:
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    //! export the scalar type
    using Scalar = ScalarType;

    //! default constructor
    CompositionalFluidState() = default;

    //! copy constructor from arbitrary fluid state
    template <class FluidState, typename std::enable_if_t<!std::is_same<FluidState, CompositionalFluidState>::value, int> = 0>
    CompositionalFluidState(const FluidState &fs)
    { assign(fs); }

    // copy and move constructor / assignment operator
    CompositionalFluidState(const CompositionalFluidState &fs) = default;
    CompositionalFluidState(CompositionalFluidState &&fs) = default;
    CompositionalFluidState& operator=(const CompositionalFluidState &fs) = default;
    CompositionalFluidState& operator=(CompositionalFluidState &&fs) = default;

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * \brief Returns the index of the most wetting phase in the
     *        fluid-solid configuration (for porous medium systems).
     */
    int wettingPhase() const { return wPhaseIdx_; }

    /*!
     * \brief Returns the saturation \f$S_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The saturation is defined as the pore space occupied by the fluid divided by the total pore space:
     *  \f[S_\alpha := \frac{\phi \mathcal{V}_\alpha}{\phi \mathcal{V}}\f]
     *
     * \param phaseIdx the index of the phase
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief Returns the molar fraction \f$x^\kappa_\alpha\f$ of the component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The molar fraction \f$x^\kappa_\alpha\f$ is defined as the ratio of the number of molecules
     * of component \f$\kappa\f$ and the total number of molecules of the phase \f$\alpha\f$.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return moleFraction_[phaseIdx][compIdx]; }

    /*!
     * \brief Returns the mass fraction \f$X^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The mass fraction \f$X^\kappa_\alpha\f$ is defined as the weight of all molecules of a
     * component divided by the total mass of the fluid phase. It is related with the component's mole fraction by means of the relation
     *
     * \f[X^\kappa_\alpha = x^\kappa_\alpha \frac{M^\kappa}{\overline M_\alpha}\;,\f]
     *
     * where \f$M^\kappa\f$ is the molar mass of component \f$\kappa\f$ and
     * \f$\overline M_\alpha\f$ is the mean molar mass of a molecule of phase
     * \f$\alpha\f$.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        // calculate the mass fractions:
        // for "mass" models this is just a back calculation
        return sumMoleFractions_[phaseIdx]
               * moleFraction(phaseIdx, compIdx)
               * FluidSystem::molarMass(compIdx)
               / averageMolarMass_[phaseIdx];
    }

    /*!
     * \brief Returns the phase mass fraction, i.e. phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     * \param phaseIdx the index of the phase
     */
    Scalar phaseMassFraction(int phaseIdx) const
    {
        Scalar totalMass = 0.0;
        for (int pIdx = 0; pIdx < numPhases; ++pIdx)
            totalMass += saturation(pIdx)*density(pIdx);

        return saturation(phaseIdx)*density(phaseIdx) / totalMass;
    }

    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * The average molar mass is the mean mass of a mole of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[\mathrm{ \overline M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa}\f]
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return averageMolarMass_[phaseIdx]; }

    /*!
     * \brief The molar concentration \f$c^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * This quantity is usually called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity \f$f^\kappa_\alpha\f$ of component \f$\kappa\f$
     *  in fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     *
     *  The fugacity is defined as:
     *  \f$f_\alpha^\kappa := \Phi^\kappa_\alpha x^\kappa_\alpha p_\alpha \;,\f$
     *  where \f$\Phi^\kappa_\alpha\f$ is the fugacity coefficient \cite reid1987 .
     *  The physical meaning of fugacity becomes clear from the equation:
     *       \f[f_\alpha^\kappa = p_\alpha \exp\left\{\frac{\zeta^\kappa_\alpha}{R T_\alpha} \right\} \;,\f]
     *  where \f$\zeta^\kappa_\alpha\f$ represents the \f$\kappa\f$'s chemical
     *  potential in phase \f$\alpha\f$, \f$R\f$ stands for the ideal gas constant,
     *  and \f$T_\alpha\f$ for the absolute temperature of phase \f$\alpha\f$. Assuming thermal equilibrium,
     *  there is a one-to-one mapping between a component's chemical potential
     *  \f$\zeta^\kappa_\alpha\f$ and its fugacity \f$f^\kappa_\alpha\f$. In this
     *  case chemical equilibrium can thus be expressed by:
     *     \f[f^\kappa := f^\kappa_\alpha = f^\kappa_\beta\quad\forall \alpha, \beta\f]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fugacityCoefficient(phaseIdx, compIdx)*moleFraction(phaseIdx, compIdx)*pressure(phaseIdx); }

    /*!
     * \brief The fugacity coefficient \f$\Phi^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fugacityCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief The molar volume \f$v_{mol,\alpha}\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[m^3/mol]}\f$
     *
     * This quantity is the inverse of the molar density.
     */
    Scalar molarVolume(int phaseIdx) const
    { return 1.0/molarDensity(phaseIdx); }

    /*!
     * \brief The mass density \f$\rho_\alpha\f$ of the fluid phase
     *  \f$\alpha\f$ in \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief The molar density \f$\rho_\alpha\f$ of the fluid phase
     *  \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*!
     * \brief The absolute temperature\f$T_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[K]}\f$
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_[phaseIdx]; }

    /*!
     * \brief The pressure \f$p_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief The partial pressure of a component in a phase \f$\mathrm{[Pa]}\f$
     * \todo is this needed?
     */
    Scalar partialPressure(int phaseIdx, int compIdx) const
    {
        assert(FluidSystem::isGas(phaseIdx));
        return moleFraction(phaseIdx, compIdx) * pressure(phaseIdx);
    }

    /*!
     * \brief The specific enthalpy \f$h_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     */
    Scalar enthalpy(int phaseIdx) const
    { return enthalpy_[phaseIdx]; }

    /*!
     * \brief The specific internal energy \f$u_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     *
     * The specific internal energy is defined by the relation:
     *
     * \f[u_\alpha = h_\alpha - \frac{p_\alpha}{\rho_\alpha}\f]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return enthalpy_[phaseIdx] - pressure(phaseIdx)/density(phaseIdx); }

    /*!
     * \brief The dynamic viscosity \f$\mu_\alpha\f$ of fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa s]}\f$
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }


    /*****************************************************
     * Access to fluid properties which only make sense
     * if assuming thermodynamic equilibrium
     *****************************************************/

    /*!
     * \brief The temperature within the domain \f$\mathrm{[K]}\f$
     */
    Scalar temperature() const
    { return temperature_[0]; }

    /*!
     * \brief The fugacity of a component  \f$\mathrm{[Pa]}\f$
     *
     * This assumes chemical equilibrium.
     */
    Scalar fugacity(int compIdx) const
    { return fugacity(0, compIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     *
     * \note If the other fluid state object is inconsistent with the
     *       thermodynamic equilibrium, the result of this method is
     *       undefined.
     */
    template <class FluidState>
    void assign(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            averageMolarMass_[phaseIdx] = 0;
            sumMoleFractions_[phaseIdx] = 0;
            temperature_[phaseIdx] = fs.temperature();
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                moleFraction_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
                fugacityCoefficient_[phaseIdx][compIdx] = fs.fugacityCoefficient(phaseIdx, compIdx);
                averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
                sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compIdx];
            }
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
            density_[phaseIdx] = fs.density(phaseIdx);
            molarDensity_[phaseIdx] = fs.molarDensity(phaseIdx);
            enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
            viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
        }
        wPhaseIdx_ = fs.wettingPhase();
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of all phases.
     */
    void setTemperature(Scalar value)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            temperature_[phaseIdx] = value;
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a specific phase.
     *        This is not implemented in this fluidstate.
     */
    void setTemperature(const int phaseIdx, const Scalar value)
    { temperature_[phaseIdx] = value; }

    /*!
     * \brief Set the fluid pressure of a phase  \f$\mathrm{[Pa]}\f$
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

    /*!
     * \brief Set the saturation of a phase \f$\mathrm{[-]}\f$
     */
    void setSaturation(int phaseIdx, Scalar value)
    { saturation_[phaseIdx] = value; }

    /*!
     * \brief Set the mole fraction of a component  in a phase \f$\mathrm{[-]}\f$
     *        and update the average molar mass \f$\mathrm{[kg/mol]}\f$ according
     *        to the current composition of the phase
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    {
        moleFraction_[phaseIdx][compIdx] = value;

        // re-calculate the mean molar mass
        sumMoleFractions_[phaseIdx] = 0.0;
        averageMolarMass_[phaseIdx] = 0.0;
        for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compJIdx];
            averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compJIdx]*FluidSystem::molarMass(compJIdx);
        }
    }

    /*!
     * \brief Set the mass fraction of a component in a phase \f$\mathrm{[-]}\f$
     *        and update the average molar mass \f$\mathrm{[kg/mol]}\f$ according
     *        to the current composition of the phase
     */
    void setMassFraction(int phaseIdx, int compIdx, Scalar value)
    {
        if (numComponents != 2)
            DUNE_THROW(Dune::NotImplemented, "This currently only works for 2 components.");
        else
        {
            // calculate average molar mass of the gas phase
            Scalar M1 = FluidSystem::molarMass(compIdx);
            Scalar M2 = FluidSystem::molarMass(1-compIdx);
            Scalar X2 = 1.0-value;
            Scalar avgMolarMass = M1*M2/(M2 + X2*(M1 - M2));

            moleFraction_[phaseIdx][compIdx] = value * avgMolarMass / M1;
            moleFraction_[phaseIdx][1-compIdx] = 1.0-moleFraction_[phaseIdx][compIdx];

            // re-calculate the mean molar mass
            sumMoleFractions_[phaseIdx] = 0.0;
            averageMolarMass_[phaseIdx] = 0.0;
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compJIdx];
                averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compJIdx]*FluidSystem::molarMass(compJIdx);
            }
        }
    }

    /*!
     * \brief Set the relative humidity of a component in a phase \f$\mathrm{[-]}\f$
     *        and update the average molar mass \f$\mathrm{[kg/mol]}\f$ according
     *        to the current composition of the phase
     */
    template <class FluidState>
    void setRelativeHumidity(FluidState &fluidState, int phaseIdx, int compIdx, Scalar value)
    {
        // asserts for the assumption under which setting the relative humidity is possible
        assert(phaseIdx == FluidSystem::nPhaseIdx);
        assert(compIdx == FluidSystem::wCompIdx);
        assert(numComponents == 2);
        assert(FluidSystem::isGas(phaseIdx));

        Scalar moleFraction = value * FluidSystem::vaporPressure(fluidState, FluidSystem::wCompIdx)
                              / fluidState.pressure(phaseIdx);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::wCompIdx, moleFraction);
        fluidState.setMoleFraction(phaseIdx, FluidSystem::nCompIdx, 1.0-moleFraction);
    }

    /*!
     * \brief Set the fugacity coefficient \f$\Phi^\kappa_\alpha\f$ of component \f$\kappa\f$
     * in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$
     */
    void setFugacityCoefficient(int phaseIdx, int compIdx, Scalar value)
    { fugacityCoefficient_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Set the density of a phase \f$\mathrm{[kg / m^3]}\f$
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Set the molar density of a phase \f$\mathrm{[mol / m^3]}\f$
     */
    void setMolarDensity(int phaseIdx, Scalar value)
    { molarDensity_[phaseIdx] = value; }

    /*!
     * \brief Set the specific enthalpy of a phase\f$\mathrm{[J/kg]}\f$
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    { enthalpy_[phaseIdx] = value; }

    /*!
     * \brief Set the dynamic viscosity of a phase \f$\mathrm{[Pa s]}\f$
     */
    void setViscosity(int phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }

    /*!
     * \brief Set the index of the most wetting phase
     */
    void setWettingPhase(int phaseIdx)
    { wPhaseIdx_ = phaseIdx; }

protected:
    //! zero-initialize all data members with braces syntax
    std::array<std::array<Scalar, numComponents>, numPhases> moleFraction_ = {};
    std::array<std::array<Scalar, numComponents>, numPhases> fugacityCoefficient_ = {};
    std::array<Scalar, numPhases> averageMolarMass_ = {};
    std::array<Scalar, numPhases> sumMoleFractions_ = {};
    std::array<Scalar, numPhases> pressure_ = {};
    std::array<Scalar, numPhases> saturation_ = {};
    std::array<Scalar, numPhases> density_ = {};
    std::array<Scalar, numPhases> molarDensity_ = {};
    std::array<Scalar, numPhases> enthalpy_ = {};
    std::array<Scalar, numPhases> viscosity_ = {};
    std::array<Scalar, numPhases> temperature_ = {};

    int wPhaseIdx_{0};
};

} // end namespace Dumux

#endif
