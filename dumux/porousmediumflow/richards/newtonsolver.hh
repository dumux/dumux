// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsModel
 * \brief A Richards model Newton solver.
 */

#ifndef DUMUX_RICHARDS_NEWTON_SOLVER_HH
#define DUMUX_RICHARDS_NEWTON_SOLVER_HH

#include <algorithm>

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
    using Indices = typename Assembler::GridVariables::VolumeVariables::Indices;
    enum { pressureIdx = Indices::pressureIdx };

    using typename ParentType::Backend;
    using typename ParentType::SolutionVector;
    using typename ParentType::ResidualVector;

public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

private:

    /*!
     * \brief Update the current solution of the Newton method
     *
     * \param varsCurrentIter The variables after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void choppedUpdate_(Variables &varsCurrentIter,
                        const SolutionVector &uLastIter,
                        const ResidualVector &deltaU) final
    {
        auto uCurrentIter = uLastIter;
        Backend::axpy(-1.0, deltaU, uCurrentIter);

        // do not clamp anything after 5 iterations
        if (this->numSteps_ <= 4)
        {
            // clamp saturation change to at most 20% per iteration
            const auto& gridGeometry = this->assembler().gridGeometry();
            auto fvGeometry = localView(gridGeometry);
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                fvGeometry.bindElement(element);

                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();

                    // calculate the old wetting phase saturation
                    const auto& spatialParams = this->assembler().problem().spatialParams();
                    const auto elemSol = elementSolution(element, uCurrentIter, gridGeometry);

                    const auto fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
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

        // update the variables
        this->solutionChanged_(varsCurrentIter, uCurrentIter);

        if (this->enableResidualCriterion())
            this->computeResidualReduction_(varsCurrentIter);
    }
};

} // end namespace Dumux

#endif
