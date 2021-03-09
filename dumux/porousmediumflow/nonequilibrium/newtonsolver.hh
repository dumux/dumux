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
 * \ingroup NonEquilibriumModel
 * \brief A MpNc specific newton solver.
 *
 * This solver calls the velocity averaging in the model after each iteration.
 */

#ifndef DUMUX_NONEQUILIBRIUM_NEWTON_SOLVER_HH
#define DUMUX_NONEQUILIBRIUM_NEWTON_SOLVER_HH

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
    using SolutionVector = typename Assembler::ResidualType;

public:
    using ParentType::ParentType;

    void newtonEndStep(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter) final
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        auto& gridVariables = this->assembler().gridVariables();
        // Averages the face velocities of a vertex. Implemented in the model.
        // The velocities are stored in the model.
        gridVariables.calcVelocityAverage(uCurrentIter);
    }
};

} // end namespace Dumux
#endif
