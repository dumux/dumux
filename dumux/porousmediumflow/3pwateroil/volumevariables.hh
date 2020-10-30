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
 * \ingroup ThreePWaterOilModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the three-phase, two-component model.
 */

#ifndef DUMUX_3P2CNI_VOLUME_VARIABLES_HH
#define DUMUX_3P2CNI_VOLUME_VARIABLES_HH

#include <vector>
#include <iostream>

#include <dumux/common/math.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>

#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

#include <dumux/common/optionalscalar.hh>
#include <dumux/common/exceptions.hh>

#include "primaryvariableswitch.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {

namespace Detail {
// helper struct and function detecting if the fluid matrix interaction features a adsorptionModel() function
#ifndef DOXYGEN // hide from doxygen
template <class FluidMatrixInteraction>
using AdsorptionModelDetector = decltype(std::declval<FluidMatrixInteraction>().adsorptionModel());
#endif // DOXYGEN

template<class FluidMatrixInteraction>
static constexpr bool hasAdsorptionModel()
{ return Dune::Std::is_detected<AdsorptionModelDetector, FluidMatrixInteraction>::value; }

}

/*!
 * \ingroup ThreePWaterOilModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the three-phase, two-component model.
 */
template <class Traits>
class ThreePWaterOilVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, ThreePWaterOilVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, ThreePWaterOilVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;
    using FS = typename Traits::FluidSystem;
    static constexpr int numFluidComps = ParentType::numFluidComponents();

    enum {
        numPs = ParentType::numFluidPhases(),

        wCompIdx = FS::wCompIdx,
        nCompIdx = FS::nCompIdx,

        wPhaseIdx = FS::wPhaseIdx,
        gPhaseIdx = FS::gPhaseIdx,
        nPhaseIdx = FS::nPhaseIdx,

        switch1Idx = ModelTraits::Indices::switch1Idx,
        switch2Idx = ModelTraits::Indices::switch2Idx,
        pressureIdx = ModelTraits::Indices::pressureIdx
    };

    // present phases
    enum {
        threePhases = ModelTraits::Indices::threePhases,
        wPhaseOnly  = ModelTraits::Indices::wPhaseOnly,
        gnPhaseOnly = ModelTraits::Indices::gnPhaseOnly,
        wnPhaseOnly = ModelTraits::Indices::wnPhaseOnly,
        gPhaseOnly  = ModelTraits::Indices::gPhaseOnly,
        wgPhaseOnly = ModelTraits::Indices::wgPhaseOnly
    };

    using EffDiffModel = typename Traits::EffectiveDiffusivityModel;
    using DiffusionCoefficients = typename Traits::DiffusionType::DiffusionCoefficientsContainer;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename Traits::FluidState;
    //! The type of the fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export the indices
    using Indices = typename ModelTraits::Indices;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! Export the primary variable switch
    using PrimaryVariableSwitch = ThreePWaterOilPrimaryVariableSwitch;
    //! State if only the gas phase is allowed to disappear
    static constexpr bool onlyGasPhaseCanDisappear()
    { return Traits::ModelTraits::onlyGasPhaseCanDisappear(); }

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
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
        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();

        if constexpr (!onlyGasPhaseCanDisappear())
        {
            /* first the saturations */
            if (phasePresence == threePhases)
            {
                sw_ = priVars[switch1Idx];
                sn_ = priVars[switch2Idx];
                sg_ = 1. - sw_ - sn_;
            }
            else if (phasePresence == wPhaseOnly)
            {
                sw_ = 1.;
                sn_ = 0.;
                sg_ = 0.;
            }
            else if (phasePresence == gnPhaseOnly)
            {
                sw_ = 0.;
                sn_ = priVars[switch1Idx];
                sg_ = 1. - sn_;
            }
            else if (phasePresence == wnPhaseOnly)
            {
                sn_ = priVars[switch2Idx];
                sw_ = 1. - sn_;
                sg_ = 0.;
            }
            else if (phasePresence == gPhaseOnly)
            {
                sw_ = 0.;
                sn_ = 0.;
                sg_ = 1.;
            }
            else if (phasePresence == wgPhaseOnly)
            {
                sw_ = priVars[switch1Idx];
                sn_ = 0.;
                sg_ = 1. - sw_;
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

            fluidState_.setSaturation(wPhaseIdx, sw_);
            fluidState_.setSaturation(gPhaseIdx, sg_);
            fluidState_.setSaturation(nPhaseIdx, sn_);

            // old material law interface is deprecated: Replace this by
            // const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);
            // after the release of 3.3, when the deprecated interface is no longer supported
            const auto fluidMatrixInteraction = Deprecated::makePcKrSw<3>(Scalar{}, problem.spatialParams(), element, scv, elemSol);

            // calculate capillary pressures
            const Scalar pcgw = fluidMatrixInteraction.pcgw(sw_, sn_);
            const Scalar pcnw = fluidMatrixInteraction.pcnw(sw_, sn_);
            const Scalar pcgn = fluidMatrixInteraction.pcgn(sw_, sn_);

            const Scalar pcAlpha = fluidMatrixInteraction.pcAlpha(sw_, sn_);
            const Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

            /* now the pressures */
            if (phasePresence == threePhases || phasePresence == gnPhaseOnly || phasePresence == gPhaseOnly || phasePresence == wgPhaseOnly)
            {
                 pg_ = priVars[pressureIdx];
                 pn_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
                 pw_ = pn_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;
            }
            else if (phasePresence == wPhaseOnly || phasePresence == wnPhaseOnly)
            {
                 pw_ = priVars[pressureIdx];
                 pn_ = pw_ + pcAlpha * pcnw + (1.-pcAlpha)*pcNW1;
                 pg_ = pn_ + pcAlpha * pcgn + (1.-pcAlpha)*(pcgw - pcNW1);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

            fluidState_.setPressure(wPhaseIdx, pw_);
            fluidState_.setPressure(gPhaseIdx, pg_);
            fluidState_.setPressure(nPhaseIdx, pn_);

            /* now the temperature */
            if (phasePresence == wPhaseOnly || phasePresence == wnPhaseOnly || phasePresence == gPhaseOnly)
            {
                 temp_ = priVars[switch1Idx];
            }
            else if (phasePresence == threePhases)
            {
                 // temp from inverse pwsat and pnsat which have to sum up to pg
                 Scalar temp = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx); // initial guess
                 for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                 {
                    fluidState_.setTemperature(phaseIdx, temp);
                 }
                 solidState_.setTemperature(temp);
                 Scalar defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                     - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                 using std::abs;
                 while(abs(defect) > 0.01) // simply a small number chosen ...
                 {
                     Scalar deltaT = 1.e-8 * temp;
                 for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                 {
                    fluidState_.setTemperature(phaseIdx, temp+deltaT);
                 }
                 solidState_.setTemperature(temp+deltaT);
                     Scalar fUp = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                      - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                  for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                  {
                     fluidState_.setTemperature(phaseIdx, temp-deltaT);
                  }
                  solidState_.setTemperature(temp-deltaT);
                     Scalar fDown = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                      - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                     temp = temp - defect * 2. * deltaT / (fUp - fDown);

                     for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                     {
                         fluidState_.setTemperature(phaseIdx, temp);
                     }
                     solidState_.setTemperature(temp);
                     defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                  - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);
                 }
                 temp_ = temp;
            }
            else if (phasePresence == wgPhaseOnly)
            {
                 // temp from inverse pwsat
                 temp_ = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx);
            }
            else if (phasePresence == gnPhaseOnly)
            {
                 // temp from inverse pnsat
                 temp_ = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, nCompIdx);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

            for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                fluidState_.setTemperature(phaseIdx, temp_);
            }
            solidState_.setTemperature(temp_);

            // now comes the tricky part: calculate phase composition
            if (phasePresence == threePhases) {

                // all phases are present, phase compositions are a
                // result of the the gas <-> liquid equilibrium.
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  pg_ - partPressH2O;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // Henry
                Scalar xwn = partPressNAPL / FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx);
                Scalar xww = 1.-xwn;

                // Not yet filled with real numbers for the NAPL phase
                Scalar xnw = partPressH2O / FluidSystem::henryCoefficient(fluidState_, nPhaseIdx,wCompIdx);
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else if (phasePresence == wPhaseOnly) {
                // only the water phase is present, water phase composition is
                // stored explicitly.

                // extract mole fractions in the water phase
                Scalar xwn = priVars[switch2Idx];
                Scalar xww = 1 - xwn;

                // write water mole fractions in the fluid state
                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xgn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;
                Scalar xgw = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx) / pg_;


                // note that the NAPL phase is actually not existing!
                // thus, this is used as phase switch criterion
                // maybe solubility would be better than this approach via Henry
                Scalar xnn = xwn * FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx) / (xgn * pg_);
                Scalar xnw = xgw*pg_ / FluidSystem::henryCoefficient(fluidState_, nPhaseIdx,wCompIdx);

                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else if (phasePresence == gnPhaseOnly) {

                // only gas and NAPL phases are present

                Scalar xnw = priVars[switch2Idx];
                Scalar xnn = 1.-xnw;
                Scalar xgn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;
                Scalar xgw = 1.-xgn;

                // note that the water phase is actually not present
                // the values are used as switching criteria
                Scalar xww = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx) / pg_;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);

            }
            else if (phasePresence == wnPhaseOnly) {
                // water and NAPL are present, phase compositions are a
                // mole fractions of non-existing gas phase are used as switching criteria
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // Henry
                Scalar xwn = partPressNAPL / FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx);
                Scalar xww = 1.-xwn;

                Scalar xnw = partPressH2O / FluidSystem::henryCoefficient(fluidState_, nPhaseIdx,wCompIdx);
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else if (phasePresence == gPhaseOnly) {
                // only the gas phase is present, gas phase composition is
                // stored explicitly here below.

                const Scalar xgn = priVars[switch2Idx];
                Scalar xgw = 1 - xgn;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);

                // note that the water and NAPL phase is actually not present
                // the values are used as switching criteria
                Scalar xww = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx) / pg_;
                Scalar xnn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else if (phasePresence == wgPhaseOnly) {
                // only water and gas phases are present
                const Scalar xgn = priVars[switch2Idx];
                Scalar xgw = 1 - xgn;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);


                Scalar xwn = xgn*pg_/FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx);
                Scalar xww = 1.-xwn;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

                // note that the NAPL phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;

                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else
                assert(false); // unhandled phase state
        } // end of if(!UseSimpleModel), i.e. the more complex version with six phase states

        else // use the simpler model with only two phase states
        {
            /* first the saturations */
            if (phasePresence == threePhases)
            {
                sw_ = priVars[switch1Idx];
                sn_ = priVars[switch2Idx];
                sg_ = 1. - sw_ - sn_;
            }
            else if (phasePresence == wnPhaseOnly)
            {
                sn_ = priVars[switch2Idx];
                sw_ = 1. - sn_;
                sg_ = 0.;
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

            fluidState_.setSaturation(wPhaseIdx, sw_);
            fluidState_.setSaturation(gPhaseIdx, sg_);
            fluidState_.setSaturation(nPhaseIdx, sn_);

            // old material law interface is deprecated: Replace this by
            // const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);
            // after the release of 3.3, when the deprecated interface is no longer supported
            const auto fluidMatrixInteraction = Deprecated::makePcKrSw<3>(Scalar{}, problem.spatialParams(), element, scv, elemSol);

            // calculate capillary pressures
            const Scalar pcgw = fluidMatrixInteraction.pcgw(sw_, sn_);
            const Scalar pcnw = fluidMatrixInteraction.pcnw(sw_, sn_);
            const Scalar pcgn = fluidMatrixInteraction.pcgn(sw_, sn_);

            const Scalar pcAlpha = fluidMatrixInteraction.pcAlpha(sw_, sn_);
            const Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

            /* now the pressures */
            if (phasePresence == threePhases)
            {
                 pg_ = priVars[pressureIdx];
                 pn_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
                 pw_ = pn_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;
            }
            else if (phasePresence == wnPhaseOnly)
            {
                 pw_ = priVars[pressureIdx];
                 pn_ = pw_ + pcAlpha * pcnw + (1.-pcAlpha)*pcNW1;
                 pg_ = pn_ + pcAlpha * pcgn + (1.-pcAlpha)*(pcgw - pcNW1);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

            fluidState_.setPressure(wPhaseIdx, pw_);
            fluidState_.setPressure(gPhaseIdx, pg_);
            fluidState_.setPressure(nPhaseIdx, pn_);

            /* now the temperature */
            if (phasePresence == wnPhaseOnly)
            {
                 temp_ = priVars[switch1Idx];
            }
            else if (phasePresence == threePhases)
            {
                 if(sn_<=1.e-10) // this threshold values is chosen arbitrarily as a small number
                 {
                     Scalar tempOnlyWater = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx);
                     temp_ = tempOnlyWater;
                 }
                 if(sw_<=1.e-10) // this threshold values is chosen arbitrarily as a small number
                 {
                     Scalar tempOnlyNAPL = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, nCompIdx);
                     temp_ = tempOnlyNAPL;
                 }
                 else
                 {
                     // temp from inverse pwsat and pnsat which have to sum up to pg
                     Scalar tempOnlyNAPL = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, nCompIdx);
                     Scalar tempOnlyWater = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx);
                     for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                     {
                        fluidState_.setTemperature(phaseIdx, tempOnlyWater);
                     }
                     solidState_.setTemperature(tempOnlyWater);
                     Scalar defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                         - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                     Scalar temp = tempOnlyWater; // initial guess
                     int counter = 0;
                     using std::abs;
                     while(abs(defect) > 0.01) // simply a small number chosen ...
                     {
                         Scalar deltaT = 1.e-6; // fixed number, but T should always be in the order of a few hundred Kelvin
                         for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                         {
                            fluidState_.setTemperature(phaseIdx, temp+deltaT);
                         }
                         solidState_.setTemperature(temp+deltaT);
                         Scalar fUp = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                          - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                        for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                        {
                            fluidState_.setTemperature(phaseIdx, temp-deltaT);
                        }
                        solidState_.setTemperature(temp-deltaT);
                         Scalar fDown = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                          - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                         temp = temp - defect * 2. * deltaT / (fUp - fDown);

                         for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                         {
                              fluidState_.setTemperature(phaseIdx, temp);
                         }
                         solidState_.setTemperature(temp);
                         defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                      - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);
                         counter +=1;
                         if (counter>10) break;
                     }
                     if ((sw_>1.e-10)&&(sw_<0.01))
                         temp = temp + (sw_ - 1.e-10) * (temp - tempOnlyNAPL) / (0.01 - 1.e-10);
                     if ((sn_>1.e-10)&&(sn_<0.01))
                         temp = temp + (sn_ - 1.e-10) * (temp - tempOnlyWater) / (0.01 - 1.e-10);
                     temp_ = temp;
                 }
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

            for (int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                fluidState_.setTemperature(phaseIdx, temp_);

            solidState_.setTemperature(temp_);

            // now comes the tricky part: calculate phase composition
            if (phasePresence == threePhases) {

                // all phases are present, phase compositions are a
                // result of the the gas <-> liquid equilibrium.
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  pg_ - partPressH2O;
                // regularized evaporation for small liquid phase saturations
                // avoids negative saturations of liquid phases
                if (sw_<0.02) partPressH2O *= sw_/0.02;
                if (partPressH2O < 0.) partPressH2O = 0;
                if (sn_<0.02) partPressNAPL *= sn_ / 0.02;
                if (partPressNAPL < 0.) partPressNAPL = 0;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // Immiscible liquid phases, mole fractions are just dummy values
                Scalar xwn = 0;
                Scalar xww = 1.-xwn;

                // Not yet filled with real numbers for the NAPL phase
                Scalar xnw = 0;
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else if (phasePresence == wnPhaseOnly) {
                // mole fractions of non-existing gas phase are used as switching criteria
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // immiscible liquid phases, mole fractions are just dummy values
                Scalar xwn = 0;
                Scalar xww = 1.-xwn;

                Scalar xnw = 0;
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
                Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
                Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
                Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
                fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
                fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
                fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        }

        // old material law interface is deprecated: Replace this by
        // const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw<3>(Scalar{}, problem.spatialParams(), element, scv, elemSol);

        for (int phaseIdx = 0; phaseIdx < numPs; ++phaseIdx)
        {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState_,
                                       phaseIdx);
            fluidState_.setViscosity(phaseIdx,mu);

            const Scalar kr = fluidMatrixInteraction.kr(phaseIdx,
                                 fluidState_.saturation(wPhaseIdx),
                                 fluidState_.saturation(nPhaseIdx));
            mobility_[phaseIdx] = kr / mu;
        }

        // material dependent parameters for NAPL adsorption (only if law is provided)
        if constexpr (Detail::hasAdsorptionModel<std::decay_t<decltype(fluidMatrixInteraction)>>())
            bulkDensTimesAdsorpCoeff_ = fluidMatrixInteraction.adsorptionModel().bulkDensTimesAdsorpCoeff();

        // porosity
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);

        auto getEffectiveDiffusionCoefficient = [&](int phaseIdx, int compIIdx, int compJIdx)
        {
            return EffDiffModel::effectiveDiffusionCoefficient(*this, phaseIdx, compIIdx, compJIdx);
        };

        effectiveDiffCoeff_.update(getEffectiveDiffusionCoefficient);

        // permeability
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        fluidState_.setTemperature(temp_);
        // the enthalpies (internal energies are directly calculated in the fluidstate
        for (int phaseIdx = 0; phaseIdx < numPs; ++phaseIdx)
        {
            Scalar h = FluidSystem::enthalpy(fluidState_, phaseIdx);
            fluidState_.setEnthalpy(phaseIdx, h);
        }

        EnergyVolVars::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Returns the phase state for the control-volume.
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
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
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
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    { return mobility_[phaseIdx]; }

    Scalar viscosity(const int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the permeability within the control volume.
     */
    Scalar permeability() const
    { return permeability_; }

    /*
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    {
        if (phaseIdx != nPhaseIdx)
            return FluidSystem::diffusionCoefficient(fluidState_, phaseIdx);
        else
            return 1.e-10;
    }

    /*!
     * \brief Returns the effective diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    { return effectiveDiffCoeff_(phaseIdx, compIIdx, compJIdx); }

    /*!
     * \brief Returns the adsorption information.
     */
    Scalar bulkDensTimesAdsorpCoeff() const
    {
        if (bulkDensTimesAdsorpCoeff_)
            return bulkDensTimesAdsorpCoeff_.value();
        else
            DUNE_THROW(Dune::NotImplemented, "Your spatialParams do not provide an adsorption model");
    }

     /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(int phaseIdx) const
    { return fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return fluidState_.enthalpy(phaseIdx); };

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    Scalar sw_, sg_, sn_, pg_, pw_, pn_, temp_;

    Scalar moleFrac_[numPs][numFluidComps];
    Scalar massFrac_[numPs][numFluidComps];

    Scalar permeability_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPs];  //!< Effective mobility within the control volume
    OptionalScalar<Scalar> bulkDensTimesAdsorpCoeff_; //!< the basis for calculating adsorbed NAPL

    //!< Binary diffusion coefficients of the 3 components in the phases
    DiffusionCoefficients effectiveDiffCoeff_;

};
} // end namespace Dumux

#endif
