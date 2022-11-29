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
 * \ingroup TwoPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase model.
 */

#ifndef DUMUX_PNM_MP_NC_VOLUME_VARIABLES_HH
#define DUMUX_PNM_MP_NC_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/mpnc/volumevariables.hh>
#include <dumux/porenetwork/mpnc/primaryvariableswitch.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup TwoPModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class Traits>
class MPNCVolumeVariables
: public Dumux::MPNCVolumeVariables<Traits>
{
    using ParentType = Dumux::MPNCVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;
    using EnergyVolVars = Dumux::EnergyVolumeVariables<Traits, Dumux::MPNCVolumeVariables<Traits>>;
    static constexpr bool enableDiffusion = ModelTraits::enableMolecularDiffusion();

    using ComponentVector = Dune::FieldVector<Scalar, ModelTraits::numFluidComponents()>;
    using CompositionFromFugacities = Dumux::CompositionFromFugacities<Scalar, typename Traits::FluidSystem>;

    using EffDiffModel = typename Traits::EffectiveDiffusivityModel;
    using DiffusionCoefficients = typename Traits::DiffusionType::DiffusionCoefficientsContainer;

public:
    //! Export type of fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of fluid state
    using FluidState = typename Traits::FluidState;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! Export the type encapsulating primary variable indices
    using Indices = typename Traits::ModelTraits::Indices;
    //! Export the primary variable switch
    using PrimaryVariableSwitch = MPNCPrimaryVariableSwitch;

    static constexpr int pwIsPrimaryVariableIndex = 0;
    static constexpr int pnIsPrimaryVariableIndex = 1;

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, ParentType::fluidState_, ParentType::solidState_);

        //calculate the remaining quantities
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(ParentType::fluidState_);

        // porosity
        updateSolidVolumeFractions(elemSol, problem, element, scv, ParentType::solidState_, ParentType::numFluidComps);

        if constexpr (enableDiffusion)
        {
            auto getEffectiveDiffusionCoefficient = [&](int phaseIdx, int compIIdx, int compJIdx)
            { return EffDiffModel::effectiveDiffusionCoefficient(*this, phaseIdx, compIIdx, compJIdx); };

            ParentType::effectiveDiffCoeff_.update(getEffectiveDiffusionCoefficient);
        }

        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, ParentType::solidState_);
        EnergyVolVars::updateEffectiveThermalConductivity();

        poreInscribedRadius_ = problem.spatialParams().poreInscribedRadius(element, scv, elemSol);
        poreVolume_ = problem.gridGeometry().poreVolume(scv.dofIndex()) * ParentType::porosity();
        surfaceTension_ = problem.spatialParams().surfaceTension(element, scv, elemSol);
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
     */

    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {

        /////////////
        // set the fluid phase temperatures
        /////////////
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        /////////////
        // set the phase saturations
        /////////////
        auto&& priVars = elemSol[scv.localDofIndex()];
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < ParentType::numFluidPhases() - 1; ++phaseIdx)
        {
            sumSat += priVars[Indices::s0Idx + phaseIdx];
            fluidState.setSaturation(phaseIdx, priVars[Indices::s0Idx + phaseIdx]);
        }
        fluidState.setSaturation(ParentType::numFluidPhases() - 1, 1.0 - sumSat);

        /////////////
        // set the phase pressures
        /////////////
        // capillary pressure parameters
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        // capillary pressures
        const auto& fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);
        const auto capPress = fluidMatrixInteraction.capillaryPressures(fluidState, wPhaseIdx);

        // add to the pressure of the first fluid phase
        const auto phaseState = priVars.state();

        // depending on which pressure is stored in the primary variables
        if(phaseState == pwIsPrimaryVariableIndex)
        {
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            const Scalar pw = priVars[Indices::p0Idx];
            for (int phaseIdx = 0; phaseIdx < ParentType::numFluidPhases(); ++phaseIdx)
                fluidState.setPressure(phaseIdx, pw - capPress[0] + capPress[phaseIdx]);
        }
        else
        {
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            const Scalar pn = priVars[Indices::p0Idx];
            for (int phaseIdx = ParentType::numFluidPhases()-1; phaseIdx >= 0; --phaseIdx)
                fluidState.setPressure(phaseIdx, pn - capPress[ParentType::numFluidPhases()-1] + capPress[phaseIdx]);
        }

        /////////////
        // set the fluid compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        ComponentVector fug;
        // retrieve component fugacities
        for (int compIdx = 0; compIdx < ParentType::numFluidComps; ++compIdx)
            fug[compIdx] = priVars[Indices::fug0Idx + compIdx];

        // calculate phase compositions
        for (int phaseIdx = 0; phaseIdx < ParentType::numFluidPhases(); ++phaseIdx) {
            // initial guess
            for (int compIdx = 0; compIdx < ParentType::numFluidComps; ++compIdx) {
                Scalar x_ij = 1.0/ParentType::numFluidComps;

                // set initial guess of the component's mole fraction
                fluidState.setMoleFraction(phaseIdx,
                                        compIdx,
                                        x_ij);
            }
            // calculate the phase composition from the component
            // fugacities
            CompositionFromFugacities::guessInitial(fluidState, paramCache, phaseIdx, fug);
            CompositionFromFugacities::solve(fluidState, paramCache, phaseIdx, fug);
        }

        // dynamic viscosities
        for (int phaseIdx = 0; phaseIdx < ParentType::numFluidPhases(); ++phaseIdx)
        {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the enthalpy
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }


    /*!
     * \brief Returns the value of the NCP-function for a phase.
     *
     *      \param phaseIdx The local index of the phases
     */
    Scalar phaseNcp(const unsigned int phaseIdx) const
    {
        //! Smooth minimum function
        auto smoothMin = [](const Scalar a, const Scalar b, const Scalar k)
        {
             using std::max; using std::min; using std::abs;
             const auto h = max(k-abs(a-b), 0.0 )/k;
             return min(a, b) - h*h*h*k*(1.0/6.0);
        };

        const Scalar aEval = ParentType::phaseNotPresentIneq(ParentType::fluidState(), phaseIdx);
        const Scalar bEval = ParentType::phasePresentIneq(ParentType::fluidState(), phaseIdx);

        const auto fisherBurmeister = [](const Scalar a, const Scalar b)
        {
            using std::sqrt;
            return sqrt(a*a + b*b) - (a + b);
        };

        static const bool useFisherBurmeister =  getParam<bool>("MPNC.FisherBurmeister", false);
        if (useFisherBurmeister)
        {
            using std::max;
            static const Scalar alpha = getParam<Scalar>("MPNC.RegularizedFisherBurmeisterAlpha", 1.0);
            return -alpha*fisherBurmeister(aEval, bEval) + (1.0 - alpha)*max(0.0, aEval)*max(0.0, bEval);
        }

        static const Scalar k = getParam<Scalar>("MPNC.SmoothMinimumK", -1);
        static const Scalar c = getParam<Scalar>("MPNC.FactorMin", 1.0);
        if (k > 0.0)
            return smoothMin(aEval, bEval, k);
        else
            return std::min(c*aEval, bEval);
    };

    Scalar poreInscribedRadius() const
    { return poreInscribedRadius_; }

    Scalar poreVolume() const
    { return poreVolume_; }

    Scalar surfaceTension() const
    { return surfaceTension_; }

    Scalar capillaryPressure() const
    {
        const auto wPhaseIdx = ParentType::fluidState().wettingPhase();
        if (ParentType::saturation(wPhaseIdx) > 1.0 - 1e-3)
            return 0.0;
        else
            return ParentType::pressure(1) - ParentType::pressure(0);
    }

protected:

    Scalar poreInscribedRadius_;
    Scalar poreVolume_;
    Scalar surfaceTension_;
};

} // end namespace Dumux

#endif
