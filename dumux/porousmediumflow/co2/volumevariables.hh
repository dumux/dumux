// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup CO2Model
 * \brief Contains the quantities which are constant within a
 *        finite volume in the CO2 model.
 */

#ifndef DUMUX_CO2_VOLUME_VARIABLES_HH
#define DUMUX_CO2_VOLUME_VARIABLES_HH

#include <array>

#include <dune/common/exceptions.hh>

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

#include "primaryvariableswitch.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup CO2Model
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the CO2 model.
 */
template <class Traits>
class TwoPTwoCCO2VolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPTwoCCO2VolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables< Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCCO2VolumeVariables<Traits> >;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;
    static constexpr int numFluidComps = ParentType::numFluidComponents();

    // component indices
    enum
    {
        comp0Idx = Traits::FluidSystem::comp0Idx,
        comp1Idx = Traits::FluidSystem::comp1Idx,
        phase0Idx = Traits::FluidSystem::phase0Idx,
        phase1Idx = Traits::FluidSystem::phase1Idx
    };

    // phase presence indices
    enum
    {
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases
    };

    // primary variable indices
    enum
    {
        switchIdx = ModelTraits::Indices::switchIdx,
        pressureIdx = ModelTraits::Indices::pressureIdx
    };

    // formulation
    static constexpr auto formulation = ModelTraits::priVarFormulation();

    // type used for the permeability
    using PermeabilityType = typename Traits::PermeabilityType;

    // type used for the diffusion coefficients
    using EffDiffModel = typename Traits::EffectiveDiffusivityModel;
    using DiffusionCoefficients = typename Traits::DiffusionType::DiffusionCoefficientsContainer;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename Traits::FluidState;
    //! The fluid system used here
    using FluidSystem = typename Traits::FluidSystem;
    //! Export the indices
    using Indices = typename ModelTraits::Indices;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! Export the type of the primary variable switch
    using PrimaryVariableSwitch = TwoPTwoCCO2PrimaryVariableSwitch;


    //! Return whether moles or masses are balanced
    static constexpr bool useMoles() { return ModelTraits::useMoles(); }
    //! Return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }

    // check for permissive combinations
    static_assert(ModelTraits::numFluidPhases() == 2, "NumPhases set in the model is not two!");
    static_assert(ModelTraits::numFluidComponents() == 2, "NumComponents set in the model is not two!");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol, const Problem& problem, const Element& element, const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        // Second instance of a parameter cache. Could be avoided if
        // diffusion coefficients also became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(), element, scv, elemSol);

        const int wPhaseIdx = fluidState_.wettingPhase();
        const int nPhaseIdx = 1 - wPhaseIdx;

        // relative permeabilities -> require wetting phase saturation as parameter!
        relativePermeability_[wPhaseIdx] = fluidMatrixInteraction.krw(saturation(wPhaseIdx));
        relativePermeability_[nPhaseIdx] = fluidMatrixInteraction.krn(saturation(wPhaseIdx));

        // porosity & permeabilty
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        auto getEffectiveDiffusionCoefficient = [&](int phaseIdx, int compIIdx, int compJIdx)
        {
            return EffDiffModel::effectiveDiffusionCoefficient(*this, phaseIdx, compIIdx, compJIdx);
        };

        effectiveDiffCoeff_.update(getEffectiveDiffusionCoefficient);

        EnergyVolVars::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Completes the fluid state.
     *
     * \note TODO: This is a lot of copy paste from the 2p2c: factor out code!
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param fluidState A container with the current (physical) state of the fluid
     * \param solidState A container with the current (physical) state of the solid
     *
     * Set temperature, saturations, capillary pressures, viscosities, densities and enthalpies.
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();

        const auto wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(), element, scv, elemSol);

        // set the saturations
        if (phasePresence == secondPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 0.0);
            fluidState.setSaturation(phase1Idx, 1.0);
        }
        else if (phasePresence == firstPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 1.0);
            fluidState.setSaturation(phase1Idx, 0.0);
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == TwoPFormulation::p0s1)
            {
                fluidState.setSaturation(phase1Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase0Idx, 1 - priVars[switchIdx]);
            }
            else
            {
                fluidState.setSaturation(phase0Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase1Idx, 1 - priVars[switchIdx]);
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase presence.");

        // set pressures of the fluid phases
        pc_ = fluidMatrixInteraction.pc(fluidState.saturation(wPhaseIdx));
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase1Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] + pc_
                                                                       : priVars[pressureIdx] - pc_);
        }
        else
        {
            fluidState.setPressure(phase1Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase0Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] - pc_
                                                                       : priVars[pressureIdx] + pc_);
        }

        // calculate the phase compositions
        typename FluidSystem::ParameterCache paramCache;
        // both phases are present
        if (phasePresence == bothPhases)
        {
            //Get the equilibrium mole fractions from the FluidSystem and set them in the fluidState
            //xCO2 = equilibrium mole fraction of CO2 in the liquid phase
            //yH2O = equilibrium mole fraction of H2O in the gas phase
            const auto xwCO2 = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, phase0Idx);
            const auto xgH2O = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, phase1Idx);
            const auto xwH2O = 1 - xwCO2;
            const auto xgCO2 = 1 - xgH2O;
            fluidState.setMoleFraction(phase0Idx, comp0Idx, xwH2O);
            fluidState.setMoleFraction(phase0Idx, comp1Idx, xwCO2);
            fluidState.setMoleFraction(phase1Idx, comp0Idx, xgH2O);
            fluidState.setMoleFraction(phase1Idx, comp1Idx, xgCO2);
        }

        // only the nonwetting phase is present, i.e. nonwetting phase
        // composition is stored explicitly.
        else if (phasePresence == secondPhaseOnly)
        {
            if( useMoles() ) // mole-fraction formulation
            {
                // set the fluid state
                fluidState.setMoleFraction(phase1Idx, comp0Idx, priVars[switchIdx]);
                fluidState.setMoleFraction(phase1Idx, comp1Idx, 1-priVars[switchIdx]);
                // TODO give values for non-existing wetting phase
                const auto xwCO2 = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, phase0Idx);
                const auto xwH2O = 1 - xwCO2;
                fluidState.setMoleFraction(phase0Idx, comp1Idx, xwCO2);
                fluidState.setMoleFraction(phase0Idx, comp0Idx, xwH2O);
            }
            else // mass-fraction formulation
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(phase1Idx, comp0Idx, priVars[switchIdx]);
                // TODO give values for non-existing wetting phase
                const auto xwCO2 = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, phase0Idx);
                const auto xwH2O = 1 - xwCO2;
                fluidState.setMoleFraction(phase0Idx, comp1Idx, xwCO2);
                fluidState.setMoleFraction(phase0Idx, comp0Idx, xwH2O);
            }
        }

        // only the wetting phase is present, i.e. wetting phase
        // composition is stored explicitly.
        else if (phasePresence == firstPhaseOnly)
        {
            if( useMoles() ) // mole-fraction formulation
            {
                // convert mass to mole fractions and set the fluid state
                fluidState.setMoleFraction(phase0Idx, comp0Idx, 1-priVars[switchIdx]);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, priVars[switchIdx]);
                //  TODO give values for non-existing nonwetting phase
                Scalar xnH2O = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, phase1Idx);
                Scalar xnCO2 = 1 - xnH2O;
                fluidState.setMoleFraction(phase1Idx, comp1Idx, xnCO2);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, xnH2O);
            }
            else // mass-fraction formulation
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(phase0Idx, comp1Idx, priVars[switchIdx]);
                //  TODO give values for non-existing nonwetting phase
                Scalar xnH2O = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, phase1Idx);
                Scalar xnCO2 = 1 - xnH2O;
                fluidState.setMoleFraction(phase1Idx, comp1Idx, xnCO2);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, xnH2O);
            }
        }

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx)
        {
            // set the viscosity and desity here if constraintsolver is not used
            paramCache.updateComposition(fluidState, phaseIdx);
            const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
            const Scalar rhoMolar = FluidSystem::molarDensity(fluidState, phaseIdx);
            fluidState.setMolarDensity(phaseIdx, rhoMolar);
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);

            // compute and set the enthalpy
            Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state within the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the average molar mass \f$\mathrm{[kg/mol]}\f$ of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.massFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mole fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume in \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(const int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[mol/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx) ; }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature within the control volume in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(const int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    { return relativePermeability_[phaseIdx]/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(phase1Idx) - fluidState_.pressure(phase0Idx); }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the average permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    {
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);
        return FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, phaseIdx, compIIdx, compJIdx);
    }

    /*!
     * \brief Returns the effective diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    { return effectiveDiffCoeff_(phaseIdx, compIIdx, compJIdx); }


    /*!
     * \brief Returns the wetting phase index
     */
    int wettingPhase() const
    { return fluidState_.wettingPhase(); }

private:
    FluidState fluidState_;
    SolidState solidState_;
    Scalar pc_;                     // The capillary pressure
    PermeabilityType permeability_; // Effective permeability within the control volume

    // Relative permeability within the control volume
    std::array<Scalar, ModelTraits::numFluidPhases()> relativePermeability_;

    // Effective diffusion coefficients for the phases
    DiffusionCoefficients effectiveDiffCoeff_;
};

} // end namespace Dumux

#endif
