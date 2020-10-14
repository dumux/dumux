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
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */

#ifndef DUMUX_2P2C_VOLUME_VARIABLES_HH
#define DUMUX_2P2C_VOLUME_VARIABLES_HH

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

#include <dumux/discretization/method.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/porousmediumflow/2pnc/primaryvariableswitch.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {

// forward declaration
template <class Traits,  bool enableChemicalNonEquilibrium, bool useConstraintSolver>
class TwoPTwoCVolumeVariablesImplementation;


/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
template <class Traits, bool useConstraintSolver = true>
using TwoPTwoCVolumeVariables =  TwoPTwoCVolumeVariablesImplementation<Traits,  Traits::ModelTraits::enableChemicalNonEquilibrium(), useConstraintSolver>;

/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 *        This is the base class for a 2p2c model with and without chemical nonequilibrium
 */
template <class Traits, class Impl>
class TwoPTwoCVolumeVariablesBase
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables< Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >;

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

    // formulations
    static constexpr auto formulation = ModelTraits::priVarFormulation();

    using PermeabilityType = typename Traits::PermeabilityType;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, typename Traits::FluidSystem>;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition< Scalar, typename Traits::FluidSystem >;
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
    //! Export the primary variable switch
    using PrimaryVariableSwitch = TwoPNCPrimaryVariableSwitch;

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
        asImp_().completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

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
     * \brief Sets complete fluid state.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
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

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(), element, scv, elemSol);

        const auto wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        // set the saturations
        if (phasePresence == firstPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 1.0);
            fluidState.setSaturation(phase1Idx, 0.0);
        }
        else if (phasePresence == secondPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 0.0);
            fluidState.setSaturation(phase1Idx, 1.0);
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
   }

    /*!
     * \brief Returns the phase state within the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control-volume.
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
    { return fluidState_.molarDensity(phaseIdx); }

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
    { return pc_; }

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

protected:
    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }
    Impl &asImp_() { return *static_cast<Impl*>(this); }
};
/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 *        Specialization for chemical equilibrium
 */
template <class Traits, bool useConstraintSolver>
class TwoPTwoCVolumeVariablesImplementation<Traits, false, useConstraintSolver>
: public TwoPTwoCVolumeVariablesBase<Traits, TwoPTwoCVolumeVariablesImplementation<Traits, false, useConstraintSolver>>
{
    using ParentType = TwoPTwoCVolumeVariablesBase< Traits, TwoPTwoCVolumeVariablesImplementation<Traits, false, useConstraintSolver>>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >;

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

    // formulations
    static constexpr auto formulation = ModelTraits::priVarFormulation();

    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, typename Traits::FluidSystem>;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition< Scalar, typename Traits::FluidSystem >;
public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename Traits::FluidState;
    //! The fluid system used here
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! Export the primary variable switch
    using PrimaryVariableSwitch = TwoPNCPrimaryVariableSwitch;

    //! Return whether moles or masses are balanced
    static constexpr bool useMoles() { return ModelTraits::useMoles(); }
    //! Return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }

    // check for permissive combinations
    static_assert(useMoles() || (!useMoles() && useConstraintSolver), "if !UseMoles, UseConstraintSolver has to be set to true");

    // The computations in the explicit composition update most probably assume a liquid-gas interface with
    // liquid as first phase. TODO: is this really needed? The constraint solver does the job anyway, doesn't it?
    static_assert(useConstraintSolver || (!FluidSystem::isGas(phase0Idx) && FluidSystem::isGas(phase1Idx)),
                   "Explicit composition calculation has to be re-checked for NON-liquid-gas equilibria");

    /*!
     * \brief Sets complete fluid state.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
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
        ParentType::completeFluidState(elemSol, problem, element, scv, fluidState, solidState);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();

        // calculate the phase compositions
        typename FluidSystem::ParameterCache paramCache;

        // If constraint solver is not used, get the phase pressures and set the fugacity coefficients here
        if(!useConstraintSolver)
        {
            for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++ phaseIdx)
            {
                assert(FluidSystem::isIdealMixture(phaseIdx));
                for (int compIdx = 0; compIdx < ModelTraits::numFluidComponents(); ++ compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }

        // now comes the tricky part: calculate phase compositions
        const Scalar p0 = fluidState.pressure(phase0Idx);
        const Scalar p1 = fluidState.pressure(phase1Idx);
        if (phasePresence == bothPhases)
        {
            // both phases are present, phase compositions are a result
            // of the equilibrium between the phases. This is the job
            // of the "MiscibleMultiPhaseComposition" constraint solver
            if(useConstraintSolver)
                MiscibleMultiPhaseComposition::solve(fluidState,
                                                     paramCache);
            // ... or calculated explicitly this way ...
            else
            {
                // get the partial pressure of the main component of the first phase within the
                // second phase == vapor pressure due to equilibrium. Note that in this case the
                // fugacityCoefficient * p is the vapor pressure (see implementation in respective fluidsystem)
                const Scalar partPressLiquid = FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx)*p0;

                // get the partial pressure of the main component of the gas phase
                const Scalar partPressGas = p1 - partPressLiquid;

                // calculate the mole fractions of the components within the nonwetting phase
                const Scalar xnn = partPressGas / p1;
                const Scalar xnw = partPressLiquid / p1;

                // calculate the mole fractions of the components within the wetting phase
                // note that in this case the fugacityCoefficient * p is the Henry Coefficient
                // (see implementation in respective fluidsystem)
                const Scalar xwn = partPressGas / (FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx)*p0);
                const Scalar xww = 1.0 - xwn;

                // set all mole fractions
                fluidState.setMoleFraction(phase0Idx, comp0Idx, xww);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, xwn);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, xnw);
                fluidState.setMoleFraction(phase1Idx, comp1Idx, xnn);
            }
        }
        else if (phasePresence == secondPhaseOnly)
        {
            // only the second phase is present, composition is stored explicitly.
            if( useMoles() )
            {
                fluidState.setMoleFraction(phase1Idx, comp1Idx, 1 - priVars[switchIdx]);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, priVars[switchIdx]);
            }
            // setMassFraction() has only to be called 1-numComponents times
            else
                fluidState.setMassFraction(phase1Idx, comp0Idx, priVars[switchIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver)
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 phase1Idx);
            // ... or calculated explicitly this way ...
            else
            {
                // note that the water phase is actually not existing!
                // thus, this is used as phase switch criterion
                const Scalar xnw = priVars[switchIdx];
                const Scalar xnn = 1.0 - xnw;

                // first, xww:
                // xnw * pn = "actual" (hypothetical) vapor pressure
                // fugacityCoefficient * pw = vapor pressure given by thermodynamic conditions
                // Here, xww is not actually the mole fraction of water in the wetting phase
                // xww is only the ratio of "actual" vapor pressure / "thermodynamic" vapor pressure
                // If xww > 1 : gas is over-saturated with water vapor,
                // condensation takes place (see switch criterion in model)
                const Scalar xww = xnw*p1/( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx)*p0 );

                // second, xwn:
                // partialPressure / xwn = Henry
                // partialPressure = xnn * pn
                // xwn = xnn * pn / Henry
                // Henry = fugacityCoefficient * pw
                const Scalar xwn = xnn*p1/( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx)*p0 );

                fluidState.setMoleFraction(phase0Idx, comp0Idx, xww);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, xwn);
            }
        }
        else if (phasePresence == firstPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            if( useMoles() ) // mole-fraction formulation
            {
                fluidState.setMoleFraction(phase0Idx, comp0Idx, 1-priVars[switchIdx]);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, priVars[switchIdx]);
            }
            // setMassFraction() has only to be called 1-numComponents times
            else // mass-fraction formulation
                fluidState.setMassFraction(phase0Idx, comp1Idx, priVars[switchIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver)
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 phase0Idx);
            // ... or calculated explicitly this way ...
            else
            {
                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                const Scalar xwn = priVars[switchIdx];

                // first, xnw:
                // psteam = xnw * pn = partial pressure of water in gas phase
                // psteam = fugacityCoefficient * pw
                const Scalar xnw = ( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx)*p0 )/p1;

                // second, xnn:
                // xwn = partialPressure / Henry
                // partialPressure = pn * xnn
                // xwn = pn * xnn / Henry
                // xnn = xwn * Henry / pn
                // Henry = fugacityCoefficient * pw
                const Scalar xnn = xwn*( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx)*p0 )/p1;

                fluidState.setMoleFraction(phase1Idx, comp1Idx, xnn);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, xnw);
            }
        }

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx)
        {
            // set the viscosity and desity here if constraintsolver is not used
            if(!useConstraintSolver)
            {
                paramCache.updateComposition(fluidState, phaseIdx);
                const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);
                Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
                fluidState.setMolarDensity(phaseIdx, rhoMolar);
            }

            // compute and set the enthalpy
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);
            Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }

};
/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 *        Specialization for chemical non-equilibrium.
 *        The equilibrium mole fraction is calculated using Henry's and Raoult's law
 */
template <class Traits, bool useConstraintSolver>
class TwoPTwoCVolumeVariablesImplementation<Traits, true, useConstraintSolver>
: public TwoPTwoCVolumeVariablesBase<Traits, TwoPTwoCVolumeVariablesImplementation<Traits, true, useConstraintSolver>>
{
    using ParentType = TwoPTwoCVolumeVariablesBase< Traits, TwoPTwoCVolumeVariablesImplementation<Traits, true, useConstraintSolver>>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >;

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

    // formulations
    static constexpr auto formulation = ModelTraits::priVarFormulation();

    using PermeabilityType = typename Traits::PermeabilityType;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, typename Traits::FluidSystem>;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition< Scalar, typename Traits::FluidSystem >;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename Traits::FluidState;
    //! The fluid system used here
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! Export the primary variable switch
    using PrimaryVariableSwitch = TwoPNCPrimaryVariableSwitch;

    //! Return whether moles or masses are balanced
    static constexpr bool useMoles() { return ModelTraits::useMoles(); }
    //! Return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }

    // check for permissive combinations
    static_assert(useMoles() || (!useMoles() && useConstraintSolver), "if !UseMoles, UseConstraintSolver has to be set to true");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    // The computations in the explicit composition update most probably assume a liquid-gas interface with
    // liquid as first phase. TODO: is this really needed? The constraint solver does the job anyway, doesn't it?
    static_assert(useConstraintSolver || (!FluidSystem::isGas(phase0Idx) && FluidSystem::isGas(phase1Idx)),
                   "Explicit composition calculation has to be re-checked for NON-liquid-gas equilibria");

    /*!
     * \brief Sets complete fluid state.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
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
        ParentType::completeFluidState(elemSol, problem, element, scv, fluidState, solidState);

        const auto& priVars = elemSol[scv.localDofIndex()];

        /////////////
        // set the fluid compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        updateMoleFraction(fluidState,
                           paramCache,
                           priVars);


        for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx)
        {
            // compute and set the enthalpy
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);
            Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }

    /*!
     * \brief Updates composition of all phases from the primary variables.
     *
     *        \param actualFluidState Container for all the secondary variables concerning the fluids
     *        \param paramCache Container for cache parameters
     *        \param priVars The primary Variables
     */
   void updateMoleFraction(FluidState &actualFluidState,
                           typename Traits::FluidSystem::ParameterCache & paramCache,
                           const typename Traits::PrimaryVariables& priVars)
   {
        const auto phasePresence = priVars.state();

        Scalar xwnNonEquil = 0.0;
        Scalar xwwNonEquil = 0.0;
        Scalar xnwNonEquil = 0.0;
        Scalar xnnNonEquil = 0.0;

        if (phasePresence == bothPhases)
        {
            xwnNonEquil = priVars[ModelTraits::numFluidPhases()];
            xwwNonEquil = 1-xwnNonEquil;
            xnwNonEquil = priVars[ModelTraits::numFluidPhases()+comp1Idx];

            //if the capillary pressure is high, we need to use Kelvin's equation to reduce the vapor pressure
            if (actualFluidState.saturation(phase0Idx) < 0.01)
            {
                const Scalar p1 = actualFluidState.pressure(phase1Idx);
                const Scalar partPressLiquid = FluidSystem::vaporPressure(actualFluidState, comp0Idx);
                xnwNonEquil =std::min(partPressLiquid/p1, xnwNonEquil);
            }
            xnnNonEquil = 1- xnwNonEquil;

            actualFluidState.setMoleFraction(phase0Idx, comp0Idx, xwwNonEquil);
            actualFluidState.setMoleFraction(phase0Idx, comp1Idx, xwnNonEquil);
            actualFluidState.setMoleFraction(phase1Idx, comp0Idx, xnwNonEquil);
            actualFluidState.setMoleFraction(phase1Idx, comp1Idx, xnnNonEquil);

            //check if we need this
            for(int phaseIdx=0; phaseIdx<ModelTraits::numFluidPhases(); ++phaseIdx)
                for (int compIdx = 0; compIdx < numFluidComps; ++compIdx)
                {
                    const Scalar phi = FluidSystem::fugacityCoefficient(actualFluidState,
                                                                        paramCache,
                                                                        phaseIdx,
                                                                        compIdx);
                    actualFluidState.setFugacityCoefficient(phaseIdx,
                                                            compIdx,
                                                            phi);
                }

            FluidState equilFluidState; // the fluidState *on the interface* i.e. chemical equilibrium
            equilFluidState.assign(actualFluidState) ;
            // If constraint solver is not used, get the phase pressures and set the fugacity coefficients here
            if(!useConstraintSolver)
            {
                for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++ phaseIdx)
                {
                    assert(FluidSystem::isIdealMixture(phaseIdx));
                    for (int compIdx = 0; compIdx < ModelTraits::numFluidComponents(); ++ compIdx) {
                        Scalar phi = FluidSystem::fugacityCoefficient(equilFluidState, paramCache, phaseIdx, compIdx);
                        equilFluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                    }
                }
            }

            // now comes the tricky part: calculate phase compositions
            const Scalar p0 = equilFluidState.pressure(phase0Idx);
            const Scalar p1 = equilFluidState.pressure(phase1Idx);
            // both phases are present, phase compositions are a result
            // of the equilibrium between the phases. This is the job
            // of the "MiscibleMultiPhaseComposition" constraint solver
            if(useConstraintSolver)
                MiscibleMultiPhaseComposition::solve(equilFluidState,
                                                     paramCache);
            // ... or calculated explicitly this way ...
            else
            {
                // get the partial pressure of the main component of the first phase within the
                // second phase == vapor pressure due to equilibrium. Note that in this case the
                // fugacityCoefficient * p is the vapor pressure (see implementation in respective fluidsystem)
                const Scalar partPressLiquid = FluidSystem::fugacityCoefficient(equilFluidState, phase0Idx, comp0Idx)*p0;

                // get the partial pressure of the main component of the gas phase
                const Scalar partPressGas = p1 - partPressLiquid;

                // calculate the mole fractions of the components within the nonwetting phase
                const Scalar xnn = partPressGas / p1;
                const Scalar xnw = partPressLiquid / p1;

                // calculate the mole fractions of the components within the wetting phase
                // note that in this case the fugacityCoefficient * p is the Henry Coefficient
                // (see implementation in respective fluidsystem)
                const Scalar xwn = partPressGas / (FluidSystem::fugacityCoefficient(equilFluidState, phase0Idx, comp1Idx)*p0);
                const Scalar xww = 1.0 - xwn;

                // set all mole fractions
                equilFluidState.setMoleFraction(phase0Idx, comp0Idx, xww);
                equilFluidState.setMoleFraction(phase0Idx, comp1Idx, xwn);
                equilFluidState.setMoleFraction(phase1Idx, comp0Idx, xnw);
                equilFluidState.setMoleFraction(phase1Idx, comp1Idx, xnn);
            }

            // Setting the equilibrium composition (in a kinetic model not necessarily the same as the actual mole fraction)
            for(int phaseIdx=0; phaseIdx<ModelTraits::numFluidPhases(); ++phaseIdx){
                for (int compIdx=0; compIdx< numFluidComps; ++ compIdx){
                    xEquil_[phaseIdx][compIdx] = equilFluidState.moleFraction(phaseIdx, compIdx);
                }
            }
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "nonequilibrium is only possible for 2 phases present ");
        }
        // compute densities of all phases
        for(int phaseIdx=0; phaseIdx<ModelTraits::numFluidPhases(); ++phaseIdx)
        {
            const Scalar rho = FluidSystem::density(actualFluidState, paramCache, phaseIdx);
            actualFluidState.setDensity(phaseIdx, rho);
            const Scalar rhoMolar = FluidSystem::molarDensity(actualFluidState, paramCache, phaseIdx);
            actualFluidState.setMolarDensity(phaseIdx, rhoMolar);
        }
   }

    /*!
     * \brief The mole fraction we would have in the case of chemical equilibrium /
     *        on the interface.
     *
     * \param phaseIdx The index of the fluid phase
     * \param compIdx The local index of the component
     */
    const Scalar xEquil(const unsigned int phaseIdx, const unsigned int compIdx) const
    {
        return xEquil_[phaseIdx][compIdx] ;
    }

private:
    std::array<std::array<Scalar, numFluidComps>, ModelTraits::numFluidPhases()> xEquil_;
};

} // end namespace Dumux

#endif
