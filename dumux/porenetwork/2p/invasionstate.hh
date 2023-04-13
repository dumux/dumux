// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief Invasion state class for the two-phase PNM.
 */
#ifndef DUMUX_PNM_2P_INVASIONSTATE_HH
#define DUMUX_PNM_2P_INVASIONSTATE_HH

#include <vector>
#include <type_traits>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/common/labels.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief This class updates the invasion state for the two-phase PNM.
 */
template<class P>
class TwoPInvasionState
{
    using Problem = P;

    template <class T>
    using GlobalCapillaryPressureDetector = decltype(std::declval<T>().globalCapillaryPressure());

    template<class T>
    static constexpr bool hasGlobalCapillaryPressure()
    { return Dune::Std::is_detected<GlobalCapillaryPressureDetector, T>::value; }

    enum class EventType {invasion, snapOff, none};

public:

    TwoPInvasionState(const Problem& problem) : problem_(problem)
    {
        // initialize the invasion state
        invadedCurrentIteration_.resize(problem.gridGeometry().gridView().size(0));
        invadedPreviousTimeStep_.resize(problem.gridGeometry().gridView().size(0));

        for (auto&& element : elements(problem.gridGeometry().gridView()))
        {
            const auto eIdx = problem.gridGeometry().elementMapper().index(element);
            invadedCurrentIteration_[eIdx] = problem.initialInvasionState(element);
            invadedPreviousTimeStep_[eIdx] = invadedCurrentIteration_[eIdx];
        }

        numThroatsInvaded_ = std::count(invadedCurrentIteration_.begin(), invadedCurrentIteration_.end(), true);
        verbose_ = getParamFromGroup<bool>(problem.paramGroup(), "InvasionState.Verbosity", true);
        restrictToGlobalCapillaryPressure_ = getParamFromGroup<bool>(problem.paramGroup(), "InvasionState.RestrictInvasionToGlobalCapillaryPressure", false);

        if constexpr (hasGlobalCapillaryPressure<Problem>())
        {
            if (restrictToGlobalCapillaryPressure_)
                std::cout << "\n *** Invasion behavior is restricted by a global capillary pressure defined in the problem! *** \n" << std::endl;
            else
                std::cout << "\n *** WARNING: global capillary pressure defined in the problem but InvasionState.RestrictInvasionToGlobalCapillaryPressure is set to false.\n"
                          << "     Invasion behavior will NOT be restricted! ***\n" << std::endl;
        }
    }

    //! Return whether a given throat is invaded or not.
    template<class Element>
    bool invaded(const Element& element) const
    {
        const auto eIdx = problem_.gridGeometry().elementMapper().index(element);
        return invadedCurrentIteration_[eIdx];
    }

    //! Return the number of currently invaded throats
    std::size_t numThroatsInvaded() const
    { return numThroatsInvaded_; }

    //! Update the invasion state of all throats. This is done after each Newton step by a call from the Newton solver.
    template<class SolutionVector, class GridVolumeVariables, class GridFluxVariablesCache>
    bool update(const SolutionVector& sol, const GridVolumeVariables& gridVolVars, GridFluxVariablesCache& gridFluxVarsCache)
    {
        hasChangedInCurrentIteration_ = false;
        auto fvGeometry = localView(problem_.gridGeometry());
        auto elemVolVars = localView(gridVolVars);
        auto elemFluxVarsCache = localView(gridFluxVarsCache);
        for (auto&& element : elements(problem_.gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                // checks if invasion or snap-off occurred after Newton iteration step
                if (const auto invasionResult = invasionSwitch_(element, elemVolVars, elemFluxVarsCache[scvf]); invasionResult)
                {
                    hasChangedInCurrentIteration_ = true;
                    if constexpr (GridFluxVariablesCache::cachingEnabled)
                    {
                        const auto eIdx = problem_.gridGeometry().elementMapper().index(element);
                        gridFluxVarsCache.cache(eIdx, scvf.index()).update(problem_, element, fvGeometry, elemVolVars, scvf, invadedCurrentIteration_[eIdx]);
                    }
                }
            }
        }
        numThroatsInvaded_ = std::count(invadedCurrentIteration_.begin(), invadedCurrentIteration_.end(), true);
        return hasChangedInCurrentIteration_;
    }

    //! Restore the old invasion state after a Newton iteration has failed.
    void reset()
    {
        hasChangedInCurrentIteration_ = false;
        invadedCurrentIteration_ = invadedPreviousTimeStep_;
    }

    //! Return whether an invasion or snap-off occurred anywhere. Can be used, e.g., for output file writing control.
    bool hasChanged() const
    { return hasChangedComparedToPreviousTimestep_; }

    //! Return whether an invasion or snap-off occurred anywhere during the current Newton iteration.
    bool hasChangedInCurrentIteration() const
    { return hasChangedInCurrentIteration_; }

    //! This is called after the Newton method has successfully finished one time step.
    void advance()
    {
        hasChangedComparedToPreviousTimestep_ = (invadedPreviousTimeStep_ != invadedCurrentIteration_);
        invadedPreviousTimeStep_ = invadedCurrentIteration_;
    }

    template<class SolutionVector, class GridVolumeVariables, class GridFluxVariablesCache>
    void checkIfCapillaryPressureIsCloseToEntryPressure(const SolutionVector& sol,
                                                        const GridVolumeVariables& gridVolVars,
                                                        const GridFluxVariablesCache& gridFluxVarsCache) const
    {
        using Scalar = typename SolutionVector::block_type::value_type;
        static const Scalar accuracyCriterion = getParamFromGroup<Scalar>(problem_.paramGroup(), "InvasionState.AccuracyCriterion", -1.0);

        if (accuracyCriterion < 0.0)
            return;

        auto fvGeometry = localView(problem_.gridGeometry());
        auto elemVolVars = localView(gridVolVars);
        auto elemFluxVarsCache = localView(gridFluxVarsCache);
        for (auto&& element : elements(problem_.gridGeometry().gridView()))
        {
            // Only consider throats which have been invaded during the current time step
            const auto eIdx = problem_.gridGeometry().elementMapper().index(element);
            if (!invadedCurrentIteration_[eIdx] || invadedPreviousTimeStep_[eIdx] == invadedCurrentIteration_[eIdx])
                continue;

            fvGeometry.bindElement(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                // checks if pc is close enough to the entry pressure value
                const auto& fluxVarsCache = elemFluxVarsCache[scvf];

                using std::max;
                const Scalar pc = max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());

                if (pc < accuracyCriterion * fluxVarsCache.pcEntry())
                    DUNE_THROW(NumericalProblem, "At element " << eIdx << ": pc " << pc << " too far away form pcEntry " << fluxVarsCache.pcEntry());
            }
        }
    }

private:

    //! The switch for determining the invasion state of a pore throat. Called at the end of each Newton step.
    template<class Element, class ElementVolumeVariables, class FluxVariablesCache>
    auto invasionSwitch_(const Element& element,
                         const ElementVolumeVariables& elemVolVars,
                         const FluxVariablesCache& fluxVarsCache)

    {
        using Scalar = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type;
        const auto& gridGeometry = problem_.gridGeometry();
        const auto& spatialParams = problem_.spatialParams();
        const auto eIdx = gridGeometry.elementMapper().index(element);
        bool invadedBeforeSwitch = invadedCurrentIteration_[eIdx];
        bool invadedAfterSwitch = invadedBeforeSwitch;

        // Result type, containing the local scv index of the pore from which the invasion/snap-off occurred
        // Evaluates to 'false' if no invasion/snap-off occurred
        struct Result
        {
            std::uint8_t localScvIdxWithCriticalPc;
            Scalar criticalPc;
            EventType event = EventType::none;

            operator bool() const
            { return event != EventType::none; }
        };

        // Block non-wetting phase flux out of the outlet
        static const auto blockNonwettingPhase = getParamFromGroup<std::vector<int>>(problem_.paramGroup(), "InvasionState.BlockNonwettingPhaseAtThroatLabel", std::vector<int>{});
        if (!blockNonwettingPhase.empty() && std::find(blockNonwettingPhase.begin(), blockNonwettingPhase.end(), gridGeometry.throatLabel(eIdx)) != blockNonwettingPhase.end())
        {
            invadedCurrentIteration_[eIdx] = false;
            return Result{}; // nothing happened
        }

        //Determine whether throat gets invaded or snap-off occurs
        const std::array<Scalar, 2> pc = { elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure() };
        const auto pcMax = std::max_element(pc.begin(), pc.end());
        const Scalar pcEntry = fluxVarsCache.pcEntry();
        const Scalar pcSnapoff = fluxVarsCache.pcSnapoff();

        // check if there is a user-specified global capillary pressure which needs to be obeyed
        if (maybeRestrictToGlobalCapillaryPressure_(pcEntry))
        {
            if (*pcMax > pcEntry)
            {
                std::cout << "Throat " << eIdx << " would have been invaded by pc of " << *pcMax << "but a global capillary pressure restricion was set in the problem.";
                std::cout << ". pcEntry: " << spatialParams.pcEntry(element, elemVolVars) << std::endl;
            }

            invadedCurrentIteration_[eIdx] = false;
            return Result{}; //nothing happened
        }

        if (*pcMax > pcEntry)
           invadedAfterSwitch = true;
        else if (*pcMax <= pcSnapoff)
           invadedAfterSwitch = false;

        invadedCurrentIteration_[eIdx] = invadedAfterSwitch;

        if (invadedBeforeSwitch == invadedAfterSwitch)
            return Result{}; // nothing happened
        else
        {
            Result result;
            result.localScvIdxWithCriticalPc = std::distance(pc.begin(), pcMax);
            result.criticalPc = *pcMax;
            result.event = !invadedBeforeSwitch && invadedAfterSwitch ? EventType::invasion : EventType::snapOff;

            if (verbose_)
            {
                const auto wPhaseIdx = spatialParams.template wettingPhase<typename ElementVolumeVariables::VolumeVariables::FluidSystem>(element, elemVolVars);
                const std::array sw = { elemVolVars[0].saturation(wPhaseIdx), elemVolVars[1].saturation(wPhaseIdx) };
                const auto vIdx = gridGeometry.gridView().indexSet().subIndex(element, result.localScvIdxWithCriticalPc, 1);
                if (result.event == EventType::invasion)
                {
                    std::cout << "Throat " << eIdx << " was invaded from pore "  << vIdx << " :";
                    std::cout << " pc: " << *pcMax;
                    std::cout << ", pcEntry: " << spatialParams.pcEntry(element, elemVolVars);
                    std::cout << ", sw: " << sw[result.localScvIdxWithCriticalPc] << std::endl;
                }
                else
                {
                    std::cout << "Snap-off occurred at throat " << eIdx << " from pore "  << vIdx << " :";
                    std::cout << " pc: " << *pcMax;
                    std::cout << ", pcSnapoff: " << spatialParams.pcSnapoff(element, elemVolVars);
                    std::cout << ", sw: " << sw[result.localScvIdxWithCriticalPc] << std::endl;
                }
            }

            return result;
        }
    }

    //! If the user has specified a global capillary pressure, check if it is lower than the given entry capillary pressure.
    //! This may be needed to exactly reproduce pc-S curves given by static network models.
    template<class Scalar>
    bool maybeRestrictToGlobalCapillaryPressure_(const Scalar pcEntry) const
    {
        if constexpr (hasGlobalCapillaryPressure<Problem>())
            return restrictToGlobalCapillaryPressure_ && (pcEntry > problem_.globalCapillaryPressure());
        else
            return false;
    }

    std::vector<bool> invadedCurrentIteration_;
    std::vector<bool> invadedPreviousTimeStep_;
    bool hasChangedInCurrentIteration_ = false;
    bool hasChangedComparedToPreviousTimestep_ = false;
    std::size_t numThroatsInvaded_;
    bool verbose_;
    bool restrictToGlobalCapillaryPressure_;

    const Problem& problem_;
};

} // end namespace Dumux::PoreNetwork

#endif
