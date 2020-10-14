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
 * \ingroup RichardsModel
 * \brief A Richards model Newton solver.
 */

#ifndef DUMUX_RICHARDS_NEWTON_SOLVER_HH
#define DUMUX_RICHARDS_NEWTON_SOLVER_HH

#include <algorithm>
#include <dumux/common/deprecated.hh>
#include <dumux/common/properties.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {
/*!
 * \ingroup RichardsModel
 * \brief A Richards model specific Newton solver.
 *
 * This solver 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton solver.
 *
 */
template <class Assembler, class LinearSolver>
class RichardsNewtonSolver : public NewtonSolver<Assembler, LinearSolver>
{
    using Scalar = typename Assembler::Scalar;
    using ParentType = NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename Assembler::ResidualType;
    using Indices = typename Assembler::GridVariables::VolumeVariables::Indices;
    enum { pressureIdx = Indices::pressureIdx };

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

        // do not clamp anything after 5 iterations
        if (this->numSteps_ <= 4)
        {
            // clamp saturation change to at most 20% per iteration
            const auto& gridGeometry = this->assembler().gridGeometry();
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                auto fvGeometry = localView(gridGeometry);
                fvGeometry.bindElement(element);

                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();

                    // calculate the old wetting phase saturation
                    const auto& spatialParams = this->assembler().problem().spatialParams();
                    const auto elemSol = elementSolution(element, uCurrentIter, gridGeometry);

                    // old material law interface is deprecated: Replace this by
                    // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
                    // after the release of 3.1, when the deprecated interface is no longer supported
                    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, spatialParams, element, scv, elemSol);

                    const Scalar pcMin = fluidMatrixInteraction.pc(1.0);
                    const Scalar pw = uLastIter[dofIdxGlobal][pressureIdx];
                    using std::max;
                    const Scalar pn = max(this->assembler().problem().nonwettingReferencePressure(), pw + pcMin);
                    const Scalar pcOld = pn - pw;
                    const Scalar SwOld = max(0.0, fluidMatrixInteraction.sw(pcOld));

                    // convert into minimum and maximum wetting phase pressures
                    const Scalar pwMin = pn - fluidMatrixInteraction.pc(SwOld - 0.2);
                    const Scalar pwMax = pn - fluidMatrixInteraction.pc(SwOld + 0.2);

                    // clamp the result
                    using std::clamp;
                    uCurrentIter[dofIdxGlobal][pressureIdx] = clamp(uCurrentIter[dofIdxGlobal][pressureIdx], pwMin, pwMax);
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
