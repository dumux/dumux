// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NonEquilibriumModel
 * \brief A MpNc specific newton solver.
 *
 * This solver calls the velocity averaging in the model after each iteration.
 */

#ifndef DUMUX_NONEQUILIBRIUM_NEWTON_SOLVER_HH
#define DUMUX_NONEQUILIBRIUM_NEWTON_SOLVER_HH

#include <dumux/common/pdesolver.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {
/*!
 * \ingroup NonEquilibriumModel
 * \brief A nonequilibrium specific newton solver.
 *
 * This solver calls the velocity averaging in the problem after each iteration.
 */
template <class Assembler, class LinearSolver>
class NonEquilibriumNewtonSolver : public NewtonSolver<Assembler, LinearSolver>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver>;

    using typename ParentType::Backend;
    using typename ParentType::SolutionVector;
    static constexpr bool assemblerExportsVariables = Detail::PDESolver::assemblerExportsVariables<Assembler>;

public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

    void newtonEndStep(Variables &varsCurrentIter,
                       const SolutionVector &uLastIter) final
    {
        ParentType::newtonEndStep(varsCurrentIter, uLastIter);
        const auto& uCurrentIter = Backend::dofs(varsCurrentIter);

        // Averages the face velocities of a vertex. Implemented in the model.
        // The velocities are stored in the model.
        if constexpr(!assemblerExportsVariables)
            this->assembler().gridVariables().calcVelocityAverage(uCurrentIter);
        else
            varsCurrentIter.calcVelocityAverage(uCurrentIter);
    }
};

} // end namespace Dumux
#endif
