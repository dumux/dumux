// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief Invasion state of the pore throats for the two-phase PNM
 */
#ifndef DUMUX_PNM_2P_INVASIONSTATE_HH
#define DUMUX_PNM_2P_INVASIONSTATE_HH

#include <vector>
#include <string>
#include <type_traits>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/tag.hh>
#include <dumux/porenetwork/common/labels.hh>

namespace Dumux::PoreNetwork {


namespace StateSwitchMethod {

//! Switched per Newton iteration, keeping the previous time step inside the invasion state itself
struct IterationWithPrevTimeStepInformation : public Utility::Tag<IterationWithPrevTimeStepInformation>
{ static std::string name() { return "iterationWithPrevTimeStepInformation"; } };

struct Iteration : public Utility::Tag<Iteration>
{ static std::string name() { return "iteration"; } };

struct EndOfTimeStep : public Utility::Tag<EndOfTimeStep>
{ static std::string name() { return "endOfTimeStep"; } };

} // end namespace StateSwitchMethod

namespace Detail {

//! The grid element type, reached through the problem's public grid geometry accessor
template<class T>
using ProblemElement = typename std::decay_t<
    decltype(std::declval<const T&>().gridGeometry().gridView())
>::template Codim<0>::Entity;

template<class T>
using ProblemInitialInvasionState = decltype(std::declval<const T&>().initialInvasionState(
    std::declval<const ProblemElement<T>&>()
));

//! whether the problem provides the initial invasion state of a throat
template<class T>
constexpr inline bool hasInitialInvasionState()
{ return Dune::Std::is_detected<ProblemInitialInvasionState, T>::value; }

template<class T>
using ProblemGlobalCapillaryPressure = decltype(std::declval<const T&>().globalCapillaryPressure());

//! whether the problem restricts the invasion by a global capillary pressure
template<class T>
constexpr inline bool hasGlobalCapillaryPressure()
{ return Dune::Std::is_detected<ProblemGlobalCapillaryPressure, T>::value; }

} // end namespace Detail

template<class P, class SwitchMethod = StateSwitchMethod::IterationWithPrevTimeStepInformation>
class TwoPInvasionState;

/*!
 * \ingroup PNMTwoPModel
 * \brief This class updates the invasion state for the two-phase PNM.
 */
template<class P>
class TwoPInvasionState<P, StateSwitchMethod::IterationWithPrevTimeStepInformation>
{
    using Problem = P;

public:
    //! export the state switch method
    using SwitchMethod = StateSwitchMethod::IterationWithPrevTimeStepInformation;

private:

    template <class T>
    using GlobalCapillaryPressureDetector = decltype(std::declval<T>().globalCapillaryPressure());

    template<class T>
    static constexpr bool hasGlobalCapillaryPressure()
    { return Dune::Std::is_detected<GlobalCapillaryPressureDetector, T>::value; }

    enum class EventType {invasion, snapOff, none};

public:

    [[deprecated("Use StateSwitchMethod::Iteration, where the previous time step is held by the problem. Will be removed after 3.12.")]]
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


/*!
 * \ingroup PNMTwoPModel
 * \brief This class updates the invasion state for the two-phase PNM.
 */
template<class P>
class TwoPInvasionState<P, StateSwitchMethod::Iteration>
{
    using Problem = P;

    enum class EventType {invasion, snapOff, none};

public:
    //! export the state method
    using SwitchMethod = StateSwitchMethod::Iteration;

    /*!
     * \brief The constructor
     * \param problem The problem
     * \param initialInvasionState A callable returning the initial invasion state of a given throat,
     *                             i.e. bool(const Element&). Defaults to all throats being non-invaded.
     */
    template<class InitialInvasionState>
    TwoPInvasionState(const Problem& problem, const InitialInvasionState& initialInvasionState)
    : problemPtr_(&problem)
    {
        // initialize the invasion state
        invaded_.resize(problem.gridGeometry().gridView().size(0));

        for (auto&& element : elements(problem.gridGeometry().gridView()))
        {
            const auto eIdx = problem.gridGeometry().elementMapper().index(element);
            invaded_[eIdx] = initialInvasionState(element);
        }

        numThroatsInvaded_ = std::count(invaded_.begin(), invaded_.end(), true);
        verbose_ = getParamFromGroup<bool>(problem.paramGroup(), "InvasionState.Verbosity", true);
        restrictToGlobalCapillaryPressure_ = getParamFromGroup<bool>(problem.paramGroup(), "InvasionState.RestrictInvasionToGlobalCapillaryPressure", false);

        if constexpr (Detail::hasGlobalCapillaryPressure<Problem>())
        {
            if (restrictToGlobalCapillaryPressure_)
                std::cout << "\n *** Invasion behavior is restricted by a global capillary pressure defined in the problem! *** \n" << std::endl;
            else
                std::cout << "\n *** WARNING: global capillary pressure defined in the problem but InvasionState.RestrictInvasionToGlobalCapillaryPressure is set to false.\n"
                          << "     Invasion behavior will NOT be restricted! ***\n" << std::endl;
        }
    }

    /*!
     * \brief Constructor taking the initial state from the problem, if it provides one
     * \note Kept so that problems implementing initialInvasionState() keep working unchanged.
     *       Prefer passing the initial state explicitly.
     */
    TwoPInvasionState(const Problem& problem)
    : TwoPInvasionState(problem, [&problem](const auto& element)
      {
          if constexpr (Detail::hasInitialInvasionState<Problem>())
              return problem.initialInvasionState(element);
          else
              return false;
      })
    {}

    //! Return whether a given throat is invaded or not.
    template<class Element>
    bool invaded(const Element& element) const
    {
        const auto eIdx = problemPtr_->gridGeometry().elementMapper().index(element);
        return invaded_[eIdx];
    }

    //! Return the number of currently invaded throats
    std::size_t numThroatsInvaded() const
    { return numThroatsInvaded_; }

    //! Update the invasion state of all throats. This is done after each Newton step by a call from the Newton solver.
    template<class SolutionVector, class GridVolumeVariables, class GridFluxVariablesCache>
    bool update(const SolutionVector& sol, const GridVolumeVariables& gridVolVars, GridFluxVariablesCache& gridFluxVarsCache)
    {
        hasChangedInCurrentIteration_ = false;
        auto fvGeometry = localView(problemPtr_->gridGeometry());
        auto elemVolVars = localView(gridVolVars);
        auto elemFluxVarsCache = localView(gridFluxVarsCache);
        for (auto&& element : elements(problemPtr_->gridGeometry().gridView()))
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
                        const auto eIdx = problemPtr_->gridGeometry().elementMapper().index(element);
                        gridFluxVarsCache.cache(eIdx, scvf.index()).update(*problemPtr_, element, fvGeometry, elemVolVars, scvf, invaded_[eIdx]);
                    }
                }
            }
        }
        numThroatsInvaded_ = std::count(invaded_.begin(), invaded_.end(), true);
        return hasChangedInCurrentIteration_;
    }

    //! Return whether an invasion or snap-off occurred anywhere during the current Newton iteration.
    bool hasChangedInCurrentIteration() const
    { return hasChangedInCurrentIteration_; }

    //! Return the invasion state of all throats, e.g. for comparing two time levels
    const std::vector<bool>& invaded() const
    { return invaded_; }


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
        const auto& gridGeometry = problemPtr_->gridGeometry();
        const auto& spatialParams = problemPtr_->spatialParams();
        const auto eIdx = gridGeometry.elementMapper().index(element);
        bool invadedBeforeSwitch = invaded_[eIdx];
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
        static const auto blockNonwettingPhase = getParamFromGroup<std::vector<int>>(problemPtr_->paramGroup(), "InvasionState.BlockNonwettingPhaseAtThroatLabel", std::vector<int>{});
        if (!blockNonwettingPhase.empty() && std::find(blockNonwettingPhase.begin(), blockNonwettingPhase.end(), gridGeometry.throatLabel(eIdx)) != blockNonwettingPhase.end())
        {
            invaded_[eIdx] = false;
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

            invaded_[eIdx] = false;
            return Result{}; //nothing happened
        }

        if (*pcMax > pcEntry)
           invadedAfterSwitch = true;
        else if (*pcMax <= pcSnapoff)
           invadedAfterSwitch = false;

        invaded_[eIdx] = invadedAfterSwitch;

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
        if constexpr (Detail::hasGlobalCapillaryPressure<Problem>())
            return restrictToGlobalCapillaryPressure_ && (pcEntry > problemPtr_->globalCapillaryPressure());
        else
            return false;
    }

    std::vector<bool> invaded_;
    bool hasChangedInCurrentIteration_ = false;
    std::size_t numThroatsInvaded_;
    bool verbose_;
    bool restrictToGlobalCapillaryPressure_;

    const Problem* problemPtr_;
};

/*!
 * \ingroup PNMTwoPModel
 * \brief This class updates the invasion state for the two-phase PNM after each time step.
 *
 * In contrast to StateSwitchMethod::Iteration, the discrete state is not switched within the
 * Newton iteration but only once the time step has converged. No additional Newton iteration is
 * therefore enforced by a switch. Two configurations are supported:
 *
 * - If the flux variables cache provides the continuous throat state theta, see
 *   TwoPRegularizedFluxVariablesCache, theta is the primary criterion. It is interpreted as the
 *   fraction of the time step at which the invasion or the snap-off of the throat occurs, i.e.
 *   the event takes place at theta*dt, and the transmissibilities are weighted with it rather
 *   than being switched discretely. The state is switched once theta exceeds a threshold.
 * - Otherwise the switch is purely based on the capillary pressure, i.e. a discrete switch that
 *   is evaluated once per time step rather than once per Newton iteration.
 */
template<class P>
class TwoPInvasionState<P, StateSwitchMethod::EndOfTimeStep>
{
    using Problem = P;

    template<class C>
    using ThroatStateThetaDetector = decltype(std::declval<C>().theta());

    //! whether a flux variables cache provides the continuous throat state theta
    template<class C>
    static constexpr bool hasThroatStateTheta()
    { return Dune::Std::is_detected<ThroatStateThetaDetector, C>::value; }

    enum class EventType {invasion, snapOff, none};

public:
    //! export the state method
    using SwitchMethod = StateSwitchMethod::EndOfTimeStep;

    /*!
     * \brief The constructor
     * \param problem The problem
     * \param initialInvasionState A callable returning the initial invasion state of a given throat,
     *                             i.e. bool(const Element&). Defaults to all throats being non-invaded.
     */
    template<class InitialInvasionState>
    TwoPInvasionState(const Problem& problem, const InitialInvasionState& initialInvasionState)
    : problemPtr_(&problem)
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
            invaded_[eIdx] = initialInvasionState(element);
        }

        numThroatsInvaded_ = std::count(invaded_.begin(), invaded_.end(), true);
        verbose_ = getParamFromGroup<bool>(problem.paramGroup(), "InvasionState.Verbosity", true);
    }

    /*!
     * \brief Constructor taking the initial state from the problem, if it provides one
     * \note Kept so that problems implementing initialInvasionState() keep working unchanged.
     *       Prefer passing the initial state explicitly.
     */
    TwoPInvasionState(const Problem& problem)
    : TwoPInvasionState(problem, [&problem](const auto& element)
      {
          if constexpr (Detail::hasInitialInvasionState<Problem>())
              return problem.initialInvasionState(element);
          else
              return false;
      })
    {}

    //! Return whether a given throat is invaded or not.
    template<class Element>
    bool invaded(const Element& element) const
    {
        const auto eIdx = problemPtr_->gridGeometry().elementMapper().index(element);
        return invaded_[eIdx];
    }

    //! Return the number of currently invaded throats
    std::size_t numThroatsInvaded() const
    { return numThroatsInvaded_; }

    //! Update the invasion state of all throats. This is done after each time step by a call from main.
    template<class SolutionVector, class GridVolumeVariables, class GridFluxVariablesCache>
    void update(const SolutionVector& sol, const GridVolumeVariables& gridVolVars, GridFluxVariablesCache& gridFluxVarsCache)
    {
        auto fvGeometry = localView(problemPtr_->gridGeometry());
        auto elemVolVars = localView(gridVolVars);
        auto elemFluxVarsCache = localView(gridFluxVarsCache);
        for (auto&& element : elements(problemPtr_->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                // checks if invasion or snap-off occurred within time step
                if (const auto invasionResult = invasionSwitch_(element, fvGeometry, elemVolVars, elemFluxVarsCache[scvf], scvf); invasionResult)
                {
                    if constexpr (GridFluxVariablesCache::cachingEnabled)
                    {
                        const auto eIdx = problemPtr_->gridGeometry().elementMapper().index(element);
                        gridFluxVarsCache.cache(eIdx, scvf.index()).update(*problemPtr_, element, fvGeometry, elemVolVars, scvf, invaded_[eIdx]);
                    }
                }
            }
        }
        numThroatsInvaded_ = std::count(invaded_.begin(), invaded_.end(), true);
    }

    //! Return the invasion state of all throats, e.g. for comparing two time levels
    const std::vector<bool>& invaded() const
    { return invaded_; }

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
        const auto& gridGeometry = problemPtr_->gridGeometry();
        const auto& spatialParams = problemPtr_->spatialParams();
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

        const auto pcMax = std::max_element(pc.begin(), pc.end());
        const Scalar pcEntry = fluxVarsCache.pcEntry();
        const Scalar pcSnapoff = fluxVarsCache.pcSnapoff();

        using std::abs;

        // If the flux variables cache carries a regularized throat state theta, the event within the
        // time step is resolved by it and theta is the primary criterion. Otherwise the switch is
        // purely based on the capillary pressure, i.e. a discrete switch evaluated once per time step.
        if constexpr (hasThroatStateTheta<FluxVariablesCache>())
        {
            const auto theta = fluxVarsCache.theta();
            if(!fluxVarsCache.invaded())
                invadedCur = theta > invasionThetaThreshold_ || *pcMax - pcEntry >  -invasionRelativePcThreshold_*abs(pcEntry);
            else
                invadedCur = theta > snapoffThetaThreshold_ || *pcMax - pcSnapoff >  -snapoffRelativePcThreshold_*abs(pcSnapoff);
        }
        else
        {
            if(!fluxVarsCache.invaded())
                invadedCur = *pcMax - pcEntry >  -invasionRelativePcThreshold_*abs(pcEntry);
            else
                invadedCur = *pcMax - pcSnapoff >  -snapoffRelativePcThreshold_*abs(pcSnapoff);
        }

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
                    if constexpr (hasThroatStateTheta<FluxVariablesCache>())
                        std::cout << ", theta: " << fluxVarsCache.theta();
                    std::cout << std::endl;
                }
                else
                {
                    std::cout << "Snap-off occurred at throat " << eIdx << " from pore "  << vIdx << " :";
                    std::cout << " pc: " << throatPc;
                    std::cout << ", pcSnapoff: " << spatialParams.pcSnapoff(element, elemVolVars);
                    std::cout << ", sw: " << sw[result.localScvIdxWithCriticalPc];
                    if constexpr (hasThroatStateTheta<FluxVariablesCache>())
                        std::cout << ", theta: " << fluxVarsCache.theta();
                    std::cout << std::endl;
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

    const Problem* problemPtr_;
};

} // end namespace Dumux::PoreNetwork

#endif
