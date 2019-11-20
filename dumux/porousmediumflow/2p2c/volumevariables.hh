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

#include <dumux/discretization/method.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/porousmediumflow/2pnc/primaryvariableswitch.hh>

namespace Dumux {

// forward declaration
template <class Traits,  bool enableChemicalNonEquilibrium>
class TwoPTwoCVolumeVariablesImplementation;


/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
template <class Traits>
using TwoPTwoCVolumeVariables =  TwoPTwoCVolumeVariablesImplementation<Traits,  Traits::ModelTraits::enableChemicalNonEquilibrium()>;

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

        // Second instance of a parameter cache. Could be avoided if
        // diffusion coefficients also became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& matParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        const int nPhaseIdx = 1 - wPhaseIdx;

        // relative permeabilities -> require wetting phase saturation as parameter!
        relativePermeability_[wPhaseIdx] = MaterialLaw::krw(matParams, saturation(wPhaseIdx));
        relativePermeability_[nPhaseIdx] = MaterialLaw::krn(matParams, saturation(wPhaseIdx));

        // binary diffusion coefficients
        diffCoeff_[phase0Idx] = FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, phase0Idx, comp0Idx, comp1Idx);
        diffCoeff_[phase1Idx] = FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, phase1Idx, comp0Idx, comp1Idx);

        // porosity & permeabilty
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
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

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
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
        pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx));
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
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        if(phaseIdx == compIdx)
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for phaseIdx = compIdx");
        else
            return diffCoeff_[phaseIdx];
    }

private:
    FluidState fluidState_;
    SolidState solidState_;

    Scalar pc_;                     // The capillary pressure
    PermeabilityType permeability_; // Effective permeability within the control volume

    // Relative permeability within the control volume
    std::array<Scalar, ModelTraits::numFluidPhases()> relativePermeability_;

    // Binary diffusion coefficients for the phases
    std::array<Scalar, ModelTraits::numFluidPhases()> diffCoeff_;

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
template <class Traits>
class TwoPTwoCVolumeVariablesImplementation<Traits, false>
: public TwoPTwoCVolumeVariablesBase<Traits, TwoPTwoCVolumeVariablesImplementation<Traits, false>>
{
    using ParentType = TwoPTwoCVolumeVariablesBase< Traits, TwoPTwoCVolumeVariablesImplementation<Traits, false>>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

    //! Export type of compositional multi-phase and single-phase flash
    using CompositionalMultiPhaseFlash = typename Traits::CompositionalMultiPhaseFlash;
    using CompositionalSinglePhaseFlash = typename Traits::CompositionalSinglePhaseFlash;

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

        // now comes the tricky part: calculate phase compositions
        if (phasePresence == bothPhases)
        {
            // both phases are present, phase compositions are a result
            // of the equilibrium between the phases. This is the job
            // of the multi-phase flashs
            CompositionalMultiPhaseFlash::solve(fluidState, paramCache);
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

            // calculate the composition of the non-existing phases (as
            // well as the densities of all phases). This is the job
            // of the single-phase flashs
            CompositionalSinglePhaseFlash::solve(fluidState, paramCache, phase1Idx);
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

            // calculate the composition of the non-existing phases (as
            // well as the densities of all phases). This is the job
            // of the single-phase flashs
            CompositionalSinglePhaseFlash::solve(fluidState, paramCache, phase0Idx);
        }

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx)
        {
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
template <class Traits>
class TwoPTwoCVolumeVariablesImplementation<Traits, true>
: public TwoPTwoCVolumeVariablesBase<Traits, TwoPTwoCVolumeVariablesImplementation<Traits, true>>
{
    using ParentType = TwoPTwoCVolumeVariablesBase< Traits, TwoPTwoCVolumeVariablesImplementation<Traits, true>>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

    //! Export type of compositional multi-phase and single-phase flash
    using CompositionalMultiPhaseFlash = typename Traits::CompositionalMultiPhaseFlash;

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

            // now comes the tricky part: calculate phase compositions
            // both phases are present, phase compositions are a result
            // of the equilibrium between the phases. This is the job
            // of the multi-phase flashs
            CompositionalMultiPhaseFlash::solve(equilFluidState, paramCache);

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
