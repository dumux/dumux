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
 * \copydoc Dumux::PNMNewtonSolver
 */
#ifndef DUMUX_MULTIDOMAIN_PNM_NEWTON_SOLVER_HH
#define DUMUX_MULTIDOMAIN_PNM_NEWTON_SOLVER_HH

#include <type_traits>
#include <iostream>
#include <dune/common/std/type_traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include "newtonconsistencychecks.hh"


namespace Dumux {

namespace Detail {
template <typename T>
using InvasionStateDetector = decltype(std::declval<T>().invasionState());

template<class T>
static constexpr bool hasInvasionState()
{ return Dune::Std::is_detected<InvasionStateDetector, T>::value; }


// Primary template
template <class VolumeVariables, typename U = int>
struct IsMPNC : std::false_type
{};

// Specialization for MPNC
template <class VolumeVariables>
struct IsMPNC <VolumeVariables, decltype((void) VolumeVariables::Indices::s0Idx, 0)> : std::true_type
{};

// Primary template
template <class VolumeVariables, typename U = int>
struct IsNonIsothermal : std::false_type
{};

// Specialization for non-isothermal models
template <class VolumeVariables>
struct IsNonIsothermal <VolumeVariables, decltype((void) VolumeVariables::Indices::temperatureIdx, 0)> : std::true_type
{};

}
/*!
 * \ingroup PNMTwoP
 * \brief A two-phase PNM specific  newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class Assembler, class LinearSolver, class CouplingManager,
          class Reassembler = DefaultPartialReassembler,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>,
          template<class, class> class NewtonConsistencyChecks = PNMNewtonConsistencyChecks>
class MultiDomainPNMNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>
{
    using ParentType =  MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>;
    using SolutionVector = typename Assembler::ResidualType;
    using Scalar = typename Assembler::Scalar;

public:
    using ParentType::ParentType;

    /*!
     * \brief Called after each Newton update
     *
     * \param uCurrentIter The current global solution vector
     * \param uLastIter The previous global solution vector
     */
     void newtonEndStep(SolutionVector &uCurrentIter,
                        const SolutionVector &uLastIter) final
    {
        // call the method of the base class
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& gridVariables = this->assembler().gridVariables(id);
            const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

            updateInvsionState_(uCurrentIter, id, std::integral_constant<bool, hasInvasionState>());
        });
    }

    /*!
     * \brief Returns true if the current solution can be considered to
     *        be accurate enough
     */
    bool newtonConverged() const final
    {
        if (Dune::any_true(switchedInLastIteration_))
            return false;

        return ParentType::newtonConverged();
    }

    void newtonFail(SolutionVector& u) final
    {
        ParentType::newtonFail(u);

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& gridVariables = this->assembler().gridVariables(id);
            const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

            resetInvasionState_(id, std::integral_constant<bool, hasInvasionState>());
        });
    }

    /*!
     * \brief Called if the Newton method ended successfully
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed() final
    {
        ParentType::newtonSucceed();

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& gridVariables = this->assembler().gridVariables(id);
            const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

            advanceInvasionState_(id, std::integral_constant<bool, hasInvasionState>());
        });
    }

private:


    /*!
     * \brief Update the current solution of the Newton method
     *
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void choppedUpdate_(SolutionVector &uCurrentIter,
                        const SolutionVector &uLastIter,
                        const SolutionVector &deltaU) final
    {
        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;


        static const int numChoppedUpdates = getParam<int>("Newton.NumChoppedUpdates", -1);

        if (numChoppedUpdates < 0 || this->numSteps_ <= numChoppedUpdates)
        {
            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
            {
                auto& gridVariables = this->assembler().gridVariables(id);
                const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

                if constexpr (hasInvasionState)
                {
                    using VolumeVariables = typename std::decay_t<decltype(gridVariables)>::GridVolumeVariables::VolumeVariables;
                    using Indices = typename VolumeVariables::Indices;

                    auto& uCurrentIterPNM = uCurrentIter[id];
                    auto& uLastIterPNM = uLastIter[id];

                    if constexpr (Detail::IsMPNC<VolumeVariables>::value)
                    {
                        std::cout << "using MPNC chopped update" << std::endl;
                        static constexpr auto numPhases = VolumeVariables::numFluidPhases();
                        static constexpr auto numComponents = VolumeVariables::numFluidComponents();
                        for (std::size_t i = 0; i < uLastIterPNM.size(); ++i)
                        {
                            for (std::size_t phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                                saturationChop_(uCurrentIterPNM[i][Indices::s0Idx + phaseIdx],
                                                uLastIterPNM[i][Indices::s0Idx + phaseIdx]);

                            pressureChop_(uCurrentIterPNM[i][Indices::p0Idx],
                                          uLastIterPNM[i][Indices::p0Idx]);

                            for (std::size_t comp = 0; comp < numComponents; ++comp)
                                pressureChop_(uCurrentIterPNM[i][Indices::fug0Idx + comp],
                                              uLastIterPNM[i][Indices::fug0Idx + comp]);

                            if constexpr (Detail::IsNonIsothermal<VolumeVariables>::value)
                                temperatureChop_(uCurrentIterPNM[i][Indices::temperatureIdx],
                                                 uLastIterPNM[i][Indices::temperatureIdx]);
                        }
                    }
                    else
                    {
                        std::cout << "using 2pnc chopped update" << std::endl;
                        for (std::size_t i = 0; i < uLastIterPNM.size(); ++i)
                        {
                            saturationChop_(uCurrentIterPNM[i][Indices::switchIdx],
                                            uLastIterPNM[i][Indices::switchIdx]);
                            pressureChop_(uCurrentIterPNM[i][Indices::pressureIdx],
                                            uLastIterPNM[i][Indices::pressureIdx]);

                            if constexpr (Detail::IsNonIsothermal<VolumeVariables>::value)
                                temperatureChop_(uCurrentIterPNM[i][Indices::temperatureIdx],
                                                 uLastIterPNM[i][Indices::temperatureIdx]);
                        }
                    }
                }
            });
        }

        if (this->enableResidualCriterion())
            this->computeResidualReduction_(uCurrentIter);

        else
        {
            // If we get here, the convergence criterion does not require
            // additional residual evalutions. Thus, the grid variables have
            // not yet been updated to the new uCurrentIter.
            this->assembler().updateGridVariables(uCurrentIter);
        }
    }

    static void clampValue_(Scalar &val,
                            const Scalar minVal,
                            const Scalar maxVal)
    {
        using std::max;
        using std::min;
        val = max(minVal, min(val, maxVal));
    };

    static void pressureChop_(Scalar &val,
                              const Scalar oldVal)
    {
        using std::max;
        static const Scalar maxDeltaInput = getParam<Scalar>("Newton.PressureMaxDelta", -1.0);
        const Scalar maxDelta = maxDeltaInput > 0.0 ? maxDeltaInput : max(oldVal/4.0, 10e3);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        val = max(0.0, val); // don't allow negative pressures
    }

    static void saturationChop_(Scalar &val,
                                const Scalar oldVal)
    {
        static const Scalar maxDelta = getParam<Scalar>("Newton.SaturationMaxDelta", 0.25);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        static const auto range = getParam<std::array<Scalar, 2>>("Newton.SaturationClampValues", std::array<Scalar, 2>{0.0, 1.0});
        clampValue_(val, range[0], range[1]);
    }

    static void temperatureChop_(Scalar &val,
                                 const Scalar oldVal)
    {
        static const Scalar maxDelta = getParam<Scalar>("Newton.TemperatureMaxDelta", 2.0);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        static const auto range = getParam<std::array<Scalar, 2>>("Newton.TemperatureClampValues", std::array<Scalar, 2>{273.0, 600.00});
        clampValue_(val, range[0], range[1]);
    }


    template<std::size_t i>
    void updateInvsionState_(const SolutionVector& uCurrentIter, Dune::index_constant<i> id, std::true_type)
    {
        std::cout << "updating invasion state for " << i << std::endl;
        auto& gridVariables = this->assembler().gridVariables(id);
        switchedInLastIteration_[i] =  gridVariables.gridFluxVarsCache().invasionState().update(uCurrentIter[id],
                                                                                                gridVariables.curGridVolVars(),
                                                                                                gridVariables.gridFluxVarsCache());

        // If the solution is about to be accepted, check for accuracy and trigger a retry
        // with a decreased time step size if necessary.
        if (newtonConverged())
        {
            NewtonConsistencyChecks<std::decay_t<decltype(gridVariables)>, std::decay_t<decltype(uCurrentIter[id])>> checks;
            checks.performChecks(gridVariables, uCurrentIter[id], this->assembler().prevSol()[id]);
        }
    }

    template<std::size_t i>
    void resetInvasionState_(Dune::index_constant<i> id, std::true_type)
    {
        this->assembler().gridVariables(id).gridFluxVarsCache().invasionState().reset();
    }

    template<std::size_t i>
    void advanceInvasionState_(Dune::index_constant<i> id, std::true_type)
    {
        this->assembler().gridVariables(id).gridFluxVarsCache().invasionState().advance();
    }

    template<std::size_t i>
    void updateInvsionState_(const SolutionVector& uCurrentIter, Dune::index_constant<i> id, std::false_type)
    {
        switchedInLastIteration_[i] = false;
    }

    template<std::size_t i>
    void resetInvasionState_(Dune::index_constant<i> id, std::false_type)
    {}

    template<std::size_t i>
    void advanceInvasionState_(Dune::index_constant<i> id, std::false_type)
    {}

    std::array<bool, Assembler::Traits::numSubDomains> switchedInLastIteration_ = {};
};
}

#endif
