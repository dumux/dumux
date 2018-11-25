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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCVolumeVariables
 */
#ifndef DUMUX_2P1C_VOLUME_VARIABLES_HH
#define DUMUX_2P1C_VOLUME_VARIABLES_HH

#include <array>

#include <dune/common/exceptions.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

#include "primaryvariableswitch.hh"

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief The volume variables (i.e. secondary variables) for the two-phase one-component model.
 */
template <class Traits>
class TwoPOneCVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPOneCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPOneCVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using FS = typename Traits::FluidSystem;
    using Idx = typename Traits::ModelTraits::Indices;
    static constexpr int numFluidComps = ParentType::numComponents();

    // primary variable indices
    enum
    {
        numPhases = Traits::ModelTraits::numPhases(),
        switchIdx = Idx::switchIdx,
        pressureIdx = Idx::pressureIdx
    };

    // component indices
    enum
    {
        comp0Idx = FS::comp0Idx,
        liquidPhaseIdx = FS::liquidPhaseIdx,
        gasPhaseIdx = FS::gasPhaseIdx
    };

    // phase presence indices
    enum
    {
        twoPhases = Idx::twoPhases,
        liquidPhaseOnly  = Idx::liquidPhaseOnly,
        gasPhaseOnly  = Idx::gasPhaseOnly,
    };

    // formulations
    static constexpr auto formulation = Traits::ModelTraits::priVarFormulation();

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename Traits::FluidState;
    //! The type of the fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! The type of the indices
    using Indices = typename Traits::ModelTraits::Indices;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! export the primary variable switch
    using PrimaryVariableSwitch = TwoPOneCPrimaryVariableSwitch;

    //! return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }

    // check for permissive combinations
    static_assert(Traits::ModelTraits::numPhases() == 2, "NumPhases set in the model is not two!");
    static_assert(Traits::ModelTraits::numComponents() == 1, "NumComponents set in the model is not one!");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
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
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // relative permeabilities
            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting phase saturation
                // as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));
            relativePermeability_[phaseIdx] = kr;
            Valgrind::CheckDefined(relativePermeability_[phaseIdx]);
        }

        // porosity & permeability
        // porosity calculation over inert volumefraction
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    }

    //! Update the fluidstate
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {

        // capillary pressure parameters
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();

        // set the saturations
        if (phasePresence == twoPhases)
        {
            if (formulation == TwoPFormulation::p0s1)
            {
                fluidState.setSaturation(gasPhaseIdx, priVars[switchIdx]);
                fluidState.setSaturation(liquidPhaseIdx, 1.0 - priVars[switchIdx]);
            }
            else
            {
                fluidState.setSaturation(liquidPhaseIdx, priVars[switchIdx]);
                fluidState.setSaturation(gasPhaseIdx, 1.0 - priVars[switchIdx]);
            }
        }
        else if (phasePresence == liquidPhaseOnly)
        {
            fluidState.setSaturation(liquidPhaseIdx, 1.0);
            fluidState.setSaturation(gasPhaseIdx, 0.0);
        }
        else if (phasePresence == gasPhaseOnly)
        {
            fluidState.setSaturation(liquidPhaseIdx, 0.0);
            fluidState.setSaturation(gasPhaseIdx, 1.0);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

        // set pressures of the fluid phases
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx));
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(liquidPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(gasPhaseIdx, (wPhaseIdx == liquidPhaseIdx) ? priVars[pressureIdx] + pc_
                                                                              : priVars[pressureIdx] - pc_);
        }
        else
        {
            fluidState.setPressure(gasPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(liquidPhaseIdx, (wPhaseIdx == liquidPhaseIdx) ? priVars[pressureIdx] - pc_
                                                                                 : priVars[pressureIdx] + pc_);
        }

        // set the temperature
        updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        // set the densities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, phaseIdx);
            Scalar rhoMolar = FluidSystem::molarDensity(fluidState, phaseIdx);

            fluidState.setDensity(phaseIdx, rho);
            fluidState.setMolarDensity(phaseIdx, rhoMolar);
        }

        //get the viscosity and mobility
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState,
                                       phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);
        }

        // the enthalpies (internal energies are directly calculated in the fluidstate
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const Scalar h = FluidSystem::enthalpy(fluidState, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    //! Depending on the phase state, the fluid temperature is either obtained as a primary variable from the solution vector
    //! or calculated from the liquid's vapor pressure.
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateTemperature(const ElemSol& elemSol,
                           const Problem& problem,
                           const Element& element,
                           const Scv& scv,
                           FluidState& fluidState,
                           SolidState& solidState)
    {
        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);

        // get temperature
        Scalar fluidTemperature;
        if (phasePresence == liquidPhaseOnly || phasePresence == gasPhaseOnly)
            fluidTemperature = priVars[switchIdx];
        else if (phasePresence == twoPhases)
            fluidTemperature = FluidSystem::vaporTemperature(fluidState, wPhaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

        Valgrind::CheckDefined(fluidTemperature);

        // the model assumes that all fluid phases have the same temperature
        for (int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            fluidState.setTemperature(phaseIdx, fluidTemperature);

        // the solid phase could have a different temperature
        if (Traits::ModelTraits::numEnergyEq() == 1)
            solidState.setTemperature(fluidTemperature);
        else
        {
            const Scalar solidTemperature = elemSol[scv.localDofIndex()][Traits::ModelTraits::numEq()-1];
            solidState.setTemperature(solidTemperature);
        }
    }

    /*!
     * \brief Returns the fluid state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }


    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the molar density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature(const int phaseIdx = 0) const
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return relativePermeability_[phaseIdx]/fluidState_.viscosity(phaseIdx);
    }

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
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the average permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the vapor temperature \f$T_{vap}(p_n)\f$ of the fluid within the control volume.
     */
    Scalar vaporTemperature() const
    { return FluidSystem::vaporTemperature(fluidState_, liquidPhaseIdx);}

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    Scalar pc_;                     //!< The capillary pressure
    PermeabilityType permeability_; //!< Effective permeability within the control volume

    //!< Relative permeability within the control volume
    std::array<Scalar, numPhases> relativePermeability_;
};

} // end namespace Dumux

#endif
