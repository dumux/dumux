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
 * \ingroup PhasefieldModel
 * \brief A phasefield model Newton solver.
 */

#ifndef DUMUX_PHASEFIELD_MULTIDOMAIN_NEWTON_SOLVER_HH
#define DUMUX_PHASEFIELD_MULTIDOMAIN_NEWTON_SOLVER_HH

#include <algorithm>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {
/*!
 * \ingroup PhasefieldModel
 * \brief A phasefield model specific Newton solver.
 *
 * This solver 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton solver.
 *
 */
template <class Assembler, class LinearSolver, class CouplingManager>
class PhasefieldMultiDomainNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver,
    CouplingManager>
{
    using Scalar = typename Assembler::Scalar;
    using ParentType = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    using SolutionVector = typename Assembler::ResidualType;

    using Indices = typename Assembler::template GridVariables<1>::VolumeVariables::Indices;

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
        auto& uCurrent = std::get<1>(uCurrentIter);

        // do not clamp anything after 5 iterations
        if (this->numSteps_ <= 4)
        {
            // clamp phasefield values to [0, 1]
            const auto domainID = Dune::index_constant<1>();
            const auto& gridGeometry = this->assembler().gridGeometry(domainID);
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                auto fvGeometry = localView(gridGeometry);
                fvGeometry.bindElement(element);

                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();

                    // calculate the old wetting phase saturation
                    //const auto& spatialParams = this->assembler().problem().spatialParams();
                    //const auto elemSol = elementSolution(element, uCurrentIter, gridGeometry);
                    //const Scalar pw = uLastIter[dofIdxGlobal][pressureIdx];
                    //using std::max;
                    //const Scalar pn = max(this->assembler().problem().nonWettingReferencePressure(), pw + pcMin);
                    //const Scalar pcOld = pn - pw;
                    //const Scalar SwOld = max(0.0, MaterialLaw::sw(materialLawParams, pcOld));

                    //// convert into minimum and maximum wetting phase pressures
                    //const Scalar pwMin = pn - MaterialLaw::pc(materialLawParams, SwOld - 0.2);
                    //const Scalar pwMax = pn - MaterialLaw::pc(materialLawParams, SwOld + 0.2);

                    // clamp the result
                    using std::clamp;
                    uCurrent[dofIdxGlobal][Indices::phi1Idx] =
                        clamp(uCurrent[dofIdxGlobal][Indices::phi1Idx], 0.0, 1.0);
                    uCurrent[dofIdxGlobal][Indices::phi2Idx] =
                        clamp(uCurrent[dofIdxGlobal][Indices::phi2Idx], 0.0, 1.0);
                    uCurrent[dofIdxGlobal][Indices::phi3Idx] =
                        clamp(uCurrent[dofIdxGlobal][Indices::phi3Idx], 0.0, 1.0);
                }
            }
        }

        // update the grid variables
        this->solutionChanged_(uCurrentIter);

        if (this->enableResidualCriterion())
            this->computeResidualReduction_(uCurrentIter);
    }
};

} // end namespace Dumux

#endif
