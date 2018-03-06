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
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

#include "indices.hh" // for formulation

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, n-component model.
 */
template <class TypeTag>
class TwoPNCVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents),

        // formulations
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        pwsn = TwoPNCFormulation::pwsn,
        pnsw = TwoPNCFormulation::pnsw,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // component indices
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        // phase presence enums
        nPhaseOnly = Indices::nPhaseOnly,
        wPhaseOnly = Indices::wPhaseOnly,
        bothPhases = Indices::bothPhases,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FluidSystem>;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static_assert(useMoles, "use moles has to be set true in the 2pnc model");

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        const auto &materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

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
            else // ATTENTION: krn requires the wetting-phase saturation as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));

            mobility_[phaseIdx] = kr / fluidState_.viscosity(phaseIdx);

            int compIIdx = phaseIdx;
            for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                // binary diffusion coefficients
                if(compIIdx!= compJIdx)
                {
                    setDiffusionCoefficient_(phaseIdx, compJIdx,
                                             FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                     paramCache,
                                                                                     phaseIdx,
                                                                                     compIIdx,
                                                                                     compJIdx));
                }
            }
        }

        // porosity & permeability
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
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
    template<class ElementSolution>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        auto&& priVars = ParentType::extractDofPriVars(elemSol, scv);
        const auto phasePresence = priVars.state();

        /////////////
        // set the saturations
        /////////////

        Scalar Sn;
        if (phasePresence == nPhaseOnly)
        {
            Sn = 1.0;
        }
        else if (phasePresence == wPhaseOnly)
        {
            Sn = 0.0;
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == pwsn)
                Sn = priVars[switchIdx];
            else if (formulation == pnsw)
                Sn = 1.0 - priVars[switchIdx];
            else
                DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        }

        fluidState.setSaturation(nPhaseIdx, Sn);
        fluidState.setSaturation(wPhaseIdx, 1.0 - Sn);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        auto pc = MaterialLaw::pc(materialParams, 1 - Sn);

        // extract the pressures
        if (formulation == pwsn)
        {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            if (priVars[pressureIdx] + pc < 0.0)
                 DUNE_THROW(NumericalProblem,"Capillary pressure is too low");
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pnsw)
        {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            // Here we check for (p_g - pc) in order to ensure that (p_l > 0)
            if (priVars[pressureIdx] - pc < 0.0)
            {
                std::cout<< "p_n: "<< priVars[pressureIdx]<<" Cap_press: "<< pc << std::endl;
                DUNE_THROW(NumericalProblem,"Capillary pressure is too high");
            }
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases)
        {
            // both phases are present, phase composition results from
            // the nonwetting <-> wetting equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver

            // set the known mole fractions in the fluidState so that they
            // can be used by the MiscibleMultiPhaseComposition constraint solver
            unsigned int knownPhaseIdx = nPhaseIdx;
            if (GET_PROP_VALUE(TypeTag, SetMoleFractionsForWettingPhase))
            {
                knownPhaseIdx = wPhaseIdx;
            }

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(knownPhaseIdx, compIdx, priVars[compIdx]);
            }

            MiscibleMultiPhaseComposition::solve(fluidState,
                                           paramCache,
                                           /*setViscosity=*/true,
                                           /*setEnthalpy=*/false,
                                           knownPhaseIdx);
        }
        else if (phasePresence == nPhaseOnly)
        {
            Dune::FieldVector<Scalar, numComponents> moleFrac;

            moleFrac[wCompIdx] = priVars[switchIdx];

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                moleFrac[compIdx] = priVars[compIdx];

            Scalar sumMoleFracOtherComponents = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                sumMoleFracOtherComponents += moleFrac[compIdx];

            sumMoleFracOtherComponents += moleFrac[wCompIdx];
            moleFrac[nCompIdx] = 1 - sumMoleFracOtherComponents;

            // Set fluid state mole fractions
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
                fluidState.setMoleFraction(nPhaseIdx, compIdx, moleFrac[compIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             nPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        else if (phasePresence == wPhaseOnly)
        {
        // only the wetting phase is present, i.e. wetting phase
        // composition is stored explicitly.
        // extract _mass_ fractions in the non-wetting phase
            Dune::FieldVector<Scalar, numComponents> moleFrac;

            moleFrac[nCompIdx] = priVars[switchIdx];

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                moleFrac[compIdx] = priVars[compIdx];

            Scalar sumMoleFracOtherComponents = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                sumMoleFracOtherComponents += moleFrac[compIdx];

            sumMoleFracOtherComponents += moleFrac[nCompIdx];
            moleFrac[wCompIdx] = 1 - sumMoleFracOtherComponents;

            // convert mass to mole fractions and set the fluid state
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, moleFrac[compIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             wPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);

            fluidState.setDensity(phaseIdx, rho);
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
    {
        if (phaseIdx < numPhases)
            return fluidState_.density(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the kinematic viscosity of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    {
        if (phaseIdx < numPhases)
            return fluidState_.viscosity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx < numPhases)
            return fluidState_.molarDensity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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
    { return fluidState_.pressure(FluidSystem::nPhaseIdx) - fluidState_.pressure(FluidSystem::wPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

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
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient called for phaseIdx = compIdx");
    }

    /*!
     * \brief Returns the molarity of a component in the phase
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     */
     Scalar molarity(int phaseIdx, int compIdx) const // [moles/m^3]
    { return fluidState_.molarity(phaseIdx, compIdx);}

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
    Scalar porosity_;               //!< Effective porosity within the control volume
    PermeabilityType permeability_; //!> Effective permeability within the control volume
    Scalar mobility_[numPhases];    //!< Effective mobility within the control volume
    FluidState fluidState_;

private:
    void setDiffusionCoefficient_(int phaseIdx, int compIdx, Scalar d)
    {
        if (compIdx < phaseIdx)
            diffCoefficient_[phaseIdx][compIdx] = std::move(d);
        else if (compIdx > phaseIdx)
            diffCoefficient_[phaseIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient for phaseIdx = compIdx doesn't exist");
    }

    std::array<std::array<Scalar, numComponents-1>, numPhases> diffCoefficient_;
};

} // end namespace Dumux

#endif
