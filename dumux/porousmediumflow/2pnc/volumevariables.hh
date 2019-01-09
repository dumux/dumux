// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup TwoPNCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, n-component model.
 */
#ifndef DUMUX_2PNC_VOLUME_VARIABLES_HH
#define DUMUX_2PNC_VOLUME_VARIABLES_HH

#include <iostream>
#include <vector>

#include <dumux/common/math.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, n-component model.
 */
template <class Traits>
class TwoPNCVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPNCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPNCVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using FS = typename Traits::FluidSystem;
    using ModelTraits = typename Traits::ModelTraits;
    static constexpr int numFluidComps = ParentType::numComponents();
    enum
    {
        numMajorComponents = ModelTraits::numPhases(),

        // phase indices
        phase0Idx = FS::phase0Idx,
        phase1Idx = FS::phase1Idx,

        // component indices
        comp0Idx = FS::comp0Idx,
        comp1Idx = FS::comp1Idx,

        // phase presence enums
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases,

        // primary variable indices
        pressureIdx = ModelTraits::Indices::pressureIdx,
        switchIdx = ModelTraits::Indices::switchIdx
    };

    static constexpr auto formulation = ModelTraits::priVarFormulation();
    static constexpr bool setFirstPhaseMoleFractions = ModelTraits::setMoleFractionsForFirstPhase();

    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FS>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FS>;


public:
    //! export fluid state type
    using FluidState = typename Traits::FluidState;
    //! export fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    //! return whether moles or masses are balanced
    static constexpr bool useMoles() { return Traits::ModelTraits::useMoles(); }
    //! return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }

    // check for permissive specifications
    static_assert(useMoles(), "use moles has to be set true in the 2pnc model");
    static_assert(ModelTraits::numPhases() == 2, "NumPhases set in the model is not two!");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& matParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        const int nPhaseIdx = 1 - wPhaseIdx;

        // mobilities -> require wetting phase saturation as parameter!
        mobility_[wPhaseIdx] = MaterialLaw::krw(matParams, saturation(wPhaseIdx))/fluidState_.viscosity(wPhaseIdx);
        mobility_[nPhaseIdx] = MaterialLaw::krn(matParams, saturation(wPhaseIdx))/fluidState_.viscosity(nPhaseIdx);

        // binary diffusion coefficients
        for (unsigned int compJIdx = 0; compJIdx < ModelTraits::numComponents(); ++compJIdx)
        {
            if(compJIdx != comp0Idx)
                setDiffusionCoefficient_( phase0Idx, compJIdx,
                                          FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                  paramCache,
                                                                                  phase0Idx,
                                                                                  comp0Idx,
                                                                                  compJIdx) );
            if(compJIdx != comp1Idx)
                setDiffusionCoefficient_( phase1Idx, compJIdx,
                                          FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                  paramCache,
                                                                                  phase1Idx,
                                                                                  comp1Idx,
                                                                                  compJIdx) );
        }

        // calculate the remaining quantities
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    }

    /*!
     * \brief Complete the fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The problem
     * \param element The element
     * \param scv The sub control volume
     * \param fluidState The fluid state
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
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

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

        // calculate the phase compositions
        typename FluidSystem::ParameterCache paramCache;

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases)
        {
            // both phases are present, phase composition results from
            // the first <-> second phase equilibrium. This is the job
            // of the "MiscibleMultiPhaseComposition" constraint solver

            // set the known mole fractions in the fluidState so that they
            // can be used by the MiscibleMultiPhaseComposition constraint solver

            const int knownPhaseIdx = setFirstPhaseMoleFractions ? phase0Idx : phase1Idx;
            for (int compIdx = numMajorComponents; compIdx < ModelTraits::numComponents(); ++compIdx)
                fluidState.setMoleFraction(knownPhaseIdx, compIdx, priVars[compIdx]);

            MiscibleMultiPhaseComposition::solve(fluidState,
                                                 paramCache,
                                                 knownPhaseIdx);
        }
        else if (phasePresence == secondPhaseOnly)
        {

            Dune::FieldVector<Scalar, ModelTraits::numComponents()> moleFrac;

            moleFrac[comp0Idx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = moleFrac[comp0Idx];

            for (int compIdx = numMajorComponents; compIdx < ModelTraits::numComponents(); ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
                sumMoleFracOtherComponents += moleFrac[compIdx];
            }

            moleFrac[comp1Idx] = 1 - sumMoleFracOtherComponents;

            // Set fluid state mole fractions
            for (int compIdx = 0; compIdx < ModelTraits::numComponents(); ++compIdx)
                fluidState.setMoleFraction(phase1Idx, compIdx, moleFrac[compIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             phase1Idx);
        }
        else if (phasePresence == firstPhaseOnly)
        {
            // only the first phase is present, i.e. first phase composition
            // is stored explicitly. extract _mass_ fractions in the second phase
            Dune::FieldVector<Scalar, ModelTraits::numComponents()> moleFrac;

            moleFrac[comp1Idx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = moleFrac[comp1Idx];
            for (int compIdx = numMajorComponents; compIdx < ModelTraits::numComponents(); ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];

                sumMoleFracOtherComponents += moleFrac[compIdx];
            }

            moleFrac[comp0Idx] = 1 - sumMoleFracOtherComponents;

            // convert mass to mole fractions and set the fluid state
            for (int compIdx = 0; compIdx < ModelTraits::numComponents(); ++compIdx)
                fluidState.setMoleFraction(phase0Idx, compIdx, moleFrac[compIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             phase0Idx);
        }
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);

            fluidState.setDensity(phaseIdx, rho);
            fluidState.setMolarDensity(phaseIdx, rhoMolar);
            fluidState.setViscosity(phaseIdx, mu);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the kinematic viscosity of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the molar density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.molarDensity(phaseIdx);

        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    Scalar averageMolarMass(int phaseIdx) const
    {
     return fluidState_.averageMolarMass(phaseIdx);
    }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return pc_; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity();  }

    /*!
     * \brief Returns the permeability within the control volume.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }


    /*!
     * \brief Returns the diffusion coefficient
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        if (compIdx < phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for phaseIdx = compIdx");
    }

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar massFraction(int phaseIdx, int compIdx) const
     { return fluidState_.massFraction(phaseIdx, compIdx); }

     /*!
      * \brief Returns the mole fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar moleFraction(int phaseIdx, int compIdx) const
     { return fluidState_.moleFraction(phaseIdx, compIdx); }

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    void setDiffusionCoefficient_(int phaseIdx, int compIdx, Scalar d)
    {
        if (compIdx < phaseIdx)
            diffCoefficient_[phaseIdx][compIdx] = std::move(d);
        else if (compIdx > phaseIdx)
            diffCoefficient_[phaseIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient for phaseIdx = compIdx doesn't exist");
    }

    Scalar pc_;                     //!< The capillary pressure
    Scalar porosity_;               //!< Effective porosity within the control volume
    PermeabilityType permeability_; //!> Effective permeability within the control volume
    Scalar mobility_[ModelTraits::numPhases()]; //!< Effective mobility within the control volume
    std::array<std::array<Scalar, ModelTraits::numComponents()-1>, ModelTraits::numPhases()> diffCoefficient_;

};

} // end namespace Dumux

#endif
