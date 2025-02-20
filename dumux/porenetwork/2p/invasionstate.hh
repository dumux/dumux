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

enum class StateSwitchMethod
{
    iteration, theta
};

template<class P, StateSwitchMethod stateMethod = StateSwitchMethod::iteration>
class TwoPInvasionState;

/*!
 * \ingroup PNMTwoPModel
 * \brief This class updates the invasion state for the two-phase PNM.
 */
template<class P>
class TwoPInvasionState<P, StateSwitchMethod::iteration>
{
    using Problem = P;

    template <class T>
    using GlobalCapillaryPressureDetector = decltype(std::declval<T>().globalCapillaryPressure());

    template<class T>
    static constexpr bool hasGlobalCapillaryPressure()
    { return Dune::Std::is_detected<GlobalCapillaryPressureDetector, T>::value; }

    enum class EventType {invasion, snapOff, none};

public:
    //! export the state method
    static constexpr StateSwitchMethod stateMethod = StateSwitchMethod::iteration;

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
        invasionRelativePcThreshold_ = getParamFromGroup<double>(problem.paramGroup(), "InvasionState.InvasionRelativePcThreshold", 1e-6);
        snapoffRelativePcThreshold_ = getParamFromGroup<double>(problem.paramGroup(), "InvasionState.SnapoffRelativePcThreshold", 1e-6);
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

    bool updateAfterTimeStep() const
    { return false; }

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

        if (*pcMax - pcEntry >  -invasionRelativePcThreshold_*pcEntry)
           invadedAfterSwitch = true;
        else if (*pcMax - pcSnapoff < snapoffRelativePcThreshold_*pcSnapoff)
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
    double invasionRelativePcThreshold_;
    double snapoffRelativePcThreshold_;

    const Problem& problem_;
};

/*!
 * \ingroup PNMTwoPModel
 * \brief This class updates the invasion state for the two-phase PNM after each time step.
 */
template<class P>
class TwoPInvasionState<P, StateSwitchMethod::theta>
{
    using Problem = P;

    enum class EventType {invasion, snapOff, none};

public:
    //! export the state method
    static constexpr StateSwitchMethod stateMethod = StateSwitchMethod::theta;

    TwoPInvasionState(const Problem& problem) : problem_(problem)
    {
        invasionThetaThreshold_ = getParamFromGroup<double>(problem.paramGroup(), "InvasionState.InvasionThetaThreshold", 1e-10);
        snapoffThetaThreshold_ = getParamFromGroup<double>(problem.paramGroup(), "InvasionState.SnapoffThetaThreshold", 1-1e-10);
        invasionRelativePcThreshold_ = getParamFromGroup<double>(problem.paramGroup(), "InvasionState.InvasionRelativePcThreshold", 1e-6);
        snapoffRelativePcThreshold_ = getParamFromGroup<double>(problem.paramGroup(), "InvasionState.SnapoffRelativePcThreshold", 1e-6);

        // initialize the invasion state
        invaded_.resize(problem.gridGeometry().gridView().size(0));

        for (auto&& element : elements(problem.gridGeometry().gridView()))
        {
            const auto eIdx = problem.gridGeometry().elementMapper().index(element);
            invaded_[eIdx] = problem.initialInvasionState(element);
        }

        numThroatsInvaded_ = std::count(invaded_.begin(), invaded_.end(), true);
        verbose_ = getParamFromGroup<bool>(problem.paramGroup(), "InvasionState.Verbosity", true);
    }

    //! Return whether a given throat is invaded or not.
    template<class Element>
    bool invaded(const Element& element) const
    {
        const auto eIdx = problem_.gridGeometry().elementMapper().index(element);
        return invaded_[eIdx];
    }

    //! Return the number of currently invaded throats
    std::size_t numThroatsInvaded() const
    { return numThroatsInvaded_; }

    //! Update the invasion state of all throats. This is done after each time step by a call from main.
    template<class SolutionVector, class GridVolumeVariables, class GridFluxVariablesCache>
    void update(const SolutionVector& sol, const GridVolumeVariables& gridVolVars, GridFluxVariablesCache& gridFluxVarsCache)
    {
        hasChanged_ = false;
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
                // checks if invasion or snap-off occurred within time step
                if (const auto invasionResult = invasionSwitch_(element, fvGeometry, elemVolVars, elemFluxVarsCache[scvf], scvf); invasionResult)
                {
                    hasChanged_ = true;
                    if constexpr (GridFluxVariablesCache::cachingEnabled)
                    {
                        const auto eIdx = problem_.gridGeometry().elementMapper().index(element);
                        gridFluxVarsCache.cache(eIdx, scvf.index()).update(problem_, element, fvGeometry, elemVolVars, scvf, invaded_[eIdx]);
                    }
                }
            }
        }
        numThroatsInvaded_ = std::count(invaded_.begin(), invaded_.end(), true);
    }

    //! Return whether an invasion or snap-off occurred anywhere. Can be used, e.g., for output file writing control.
    bool hasChanged() const
    { return hasChanged_; }

    bool updateAfterTimeStep() const
    { return true; }

private:

    //! The switch for determining the invasion state of a pore throat. Called at the end of time step.
    template<class Element, class FvGeometry, class ElementVolumeVariables, class FluxVariablesCache, class SubControlVolumeFace>
    auto invasionSwitch_(const Element& element,
                         const FvGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const FluxVariablesCache& fluxVarsCache,
                         const SubControlVolumeFace& scvf)

    {
        using Scalar = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type;
        const auto& gridGeometry = problem_.gridGeometry();
        const auto& spatialParams = problem_.spatialParams();
        const auto eIdx = gridGeometry.elementMapper().index(element);
        bool invadedPrev = invaded_[eIdx];
        bool invadedCur = invadedPrev;

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

        // Determine whether throat gets invaded or snap-off occurs
        const std::array<Scalar, 2> pc = { elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure() };
        const auto throatPc = fluxVarsCache.pc();
        const auto theta = problem_.theta(element, fvGeometry, elemVolVars, fluxVarsCache, scvf);

        const auto pcMax = std::max_element(pc.begin(), pc.end());
        const Scalar pcEntry = fluxVarsCache.pcEntry();
        const Scalar pcSnapoff = fluxVarsCache.pcSnapoff();

        if(!fluxVarsCache.invaded())
            invadedCur = theta > invasionThetaThreshold_ || *pcMax - pcEntry >  -invasionRelativePcThreshold_*pcEntry;
        else
            invadedCur = theta > snapoffThetaThreshold_ || *pcMax - pcSnapoff <  snapoffRelativePcThreshold_*pcSnapoff;

        invaded_[eIdx] = invadedCur;

        if (invadedPrev == invadedCur)
            return Result{}; // nothing happened
        else
        {
            Result result;
            result.localScvIdxWithCriticalPc = std::abs(pc[0] - throatPc) < std::abs(pc[1] - throatPc) ? 0 : 1;
            result.criticalPc = throatPc;
            result.event = !invadedPrev && invadedCur ? EventType::invasion : EventType::snapOff;

            if (verbose_)
            {
                const auto wPhaseIdx = spatialParams.template wettingPhase<typename ElementVolumeVariables::VolumeVariables::FluidSystem>(element, elemVolVars);
                const std::array sw = { elemVolVars[0].saturation(wPhaseIdx), elemVolVars[1].saturation(wPhaseIdx) };
                const auto vIdx = gridGeometry.gridView().indexSet().subIndex(element, result.localScvIdxWithCriticalPc, 1);
                if (result.event == EventType::invasion)
                {
                    std::cout << "Throat " << eIdx << " was invaded from pore "  << vIdx << " :";
                    std::cout << " pc: " << throatPc;
                    std::cout << ", pcEntry: " << spatialParams.pcEntry(element, elemVolVars);
                    std::cout << ", sw: " << sw[result.localScvIdxWithCriticalPc];
                    std::cout << ", theta: " << theta << std::endl;
                }
                else
                {
                    std::cout << "Snap-off occurred at throat " << eIdx << " from pore "  << vIdx << " :";
                    std::cout << " pc: " << throatPc;
                    std::cout << ", pcSnapoff: " << spatialParams.pcSnapoff(element, elemVolVars);
                    std::cout << ", sw: " << sw[result.localScvIdxWithCriticalPc];
                    std::cout << ", theta: " << theta << std::endl;
                }
            }

            return result;
        }
    }

    double invasionThetaThreshold_;
    double snapoffThetaThreshold_;
    double invasionRelativePcThreshold_;
    double snapoffRelativePcThreshold_;
    std::vector<bool> invaded_;
    std::size_t numThroatsInvaded_;
    bool verbose_;
    bool hasChanged_ = false;

    const Problem& problem_;
};

} // end namespace Dumux::PoreNetwork

#endif
