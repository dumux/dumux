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
 * \copydoc Dumux::PNMNewtonSolver
 */
#ifndef DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PNM_NEWTON_SOLVER_HH
#define DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PNM_NEWTON_SOLVER_HH

#include <type_traits>
#include <iostream>
#include <dune/common/std/type_traits.hh>
#include <dumux/multidomain/newtonsolver.hh>


namespace Dumux {


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
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>>
class MultiDomainPNMNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>
{
    using ParentType =  MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>;
    using SolutionVector = typename Assembler::ResidualType;
    using Scalar = typename Assembler::Scalar;

public:
    using ParentType::ParentType;


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

                using VolumeVariables = typename std::decay_t<decltype(gridVariables)>::GridVolumeVariables::VolumeVariables;
                using Indices = typename VolumeVariables::Indices;

                auto& uCurrentIterPNM = uCurrentIter[id];
                auto& uLastIterPNM = uLastIter[id];


                std::cout << "using crazy chopped update" << std::endl;
                for (std::size_t i = 0; i < uLastIterPNM.size(); ++i)
                {
                    if constexpr (VolumeVariables::PrimaryVariables::dimension > 1)
                        pressureChop_(uCurrentIterPNM[i][Indices::pressureIdx],
                                        uLastIterPNM[i][Indices::pressureIdx]);

                    temperatureChop_(uCurrentIterPNM[i][Indices::temperatureIdx],
                                     uLastIterPNM[i][Indices::temperatureIdx]);

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


    static void temperatureChop_(Scalar &val,
                                 const Scalar oldVal)
    {
        static const Scalar maxDelta = getParam<Scalar>("Newton.TemperatureMaxDelta", 2.0);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        static const auto range = getParam<std::array<Scalar, 2>>("Newton.TemperatureClampValues", std::array<Scalar, 2>{273.0, 600.00});
        clampValue_(val, range[0], range[1]);
    }
};
}

#endif
