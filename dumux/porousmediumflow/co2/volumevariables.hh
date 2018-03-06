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
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the CO2 model.
 */
#ifndef DUMUX_CO2_VOLUME_VARIABLES_HH
#define DUMUX_CO2_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p2c/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup CO2Model
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the CO2 model.
 */
template <class TypeTag>
class TwoPTwoCCO2VolumeVariables : public TwoPTwoCVolumeVariables<TypeTag>
{
    using ParentType = TwoPTwoCVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    enum {
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    // present phases
    enum {
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    // formulations
    enum {
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        pwsn = TwoPTwoCFormulation::pwsn,
        pnsw = TwoPTwoCFormulation::pnsw
    };

    // primary variable indices
    enum {
        switchIdx = Indices::switchIdx,
        pressureIdx = Indices::pressureIdx
    };

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    //! TODO: This is a lot of copy paste from the 2p2c: factor out code!
    template<class ElementSolution>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {

        const Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        const auto phasePresence = priVars.state();

        /////////////
        // set the saturations
        /////////////
        Scalar sn;
        if (phasePresence == nPhaseOnly)
            sn = 1.0;
        else if (phasePresence == wPhaseOnly) {
            sn = 0.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == pwsn)
                sn = priVars[switchIdx];
            else if (formulation == pnsw)
                sn = 1.0 - priVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        fluidState.setSaturation(wPhaseIdx, 1 - sn);
        fluidState.setSaturation(nPhaseIdx, sn);

        // capillary pressure parameters
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        Scalar pc = MaterialLaw::pc(materialParams, 1 - sn);

        if (formulation == pwsn) {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pnsw) {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        /////////////
        // calculate the phase compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        // both phases are present
        if (phasePresence == bothPhases)
        {
            //Get the equilibrium mole fractions from the FluidSystem and set them in the fluidState
            //xCO2 = equilibrium mole fraction of CO2 in the liquid phase
            //yH2O = equilibrium mole fraction of H2O in the gas phase
            const auto xwCO2 = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, wPhaseIdx);
            const auto xgH2O = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, nPhaseIdx);
            const auto xwH2O = 1 - xwCO2;
            const auto xgCO2 = 1 - xgH2O;
            fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xwH2O);
            fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwCO2);
            fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xgH2O);
            fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xgCO2);
        }

        // only the nonwetting phase is present, i.e. nonwetting phase
        // composition is stored explicitly.
        else if (phasePresence == nPhaseOnly)
        {
            if(useMoles) // mole-fraction formulation
            {
                // set the fluid state
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, priVars[switchIdx]);
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1-priVars[switchIdx]);
                // TODO give values for non-existing wetting phase
                const auto xwCO2 = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, wPhaseIdx);
                const auto xwH2O = 1 - xwCO2;
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwCO2);
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xwH2O);
            }
            else // mass-fraction formulation
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(nPhaseIdx, wCompIdx, priVars[switchIdx]);
                // TODO give values for non-existing wetting phase
                const auto xwCO2 = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, wPhaseIdx);
                const auto xwH2O = 1 - xwCO2;
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwCO2);
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xwH2O);
            }
        }

        // only the wetting phase is present, i.e. wetting phase
        // composition is stored explicitly.
        else if (phasePresence == wPhaseOnly)
        {
            if(useMoles) // mole-fraction formulation
            {
                // convert mass to mole fractions and set the fluid state
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1-priVars[switchIdx]);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, priVars[switchIdx]);
                //  TODO give values for non-existing nonwetting phase
                Scalar xnH2O = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, nPhaseIdx);
                Scalar xnCO2 = 1 - xnH2O;
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xnCO2);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnH2O);
            }
            else // mass-fraction formulation
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(wPhaseIdx, nCompIdx, priVars[switchIdx]);
                //  TODO give values for non-existing nonwetting phase
                Scalar xnH2O = FluidSystem::equilibriumMoleFraction(fluidState, paramCache, nPhaseIdx);
                Scalar xnCO2 = 1 - xnH2O;
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xnCO2);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnH2O);
            }
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // set the viscosity and desity here if constraintsolver is not used
            paramCache.updateComposition(fluidState, phaseIdx);
            const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);

            // compute and set the enthalpy
            Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }
};

} // end namespace Dumux

#endif
