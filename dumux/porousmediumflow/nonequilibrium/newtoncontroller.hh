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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief A MpNc specific controller for the newton solver.
 *  This controller calls the velocity averaging in the model after each iteration.
 */
#ifndef DUMUX_NONEQUILIBRIUM_NEWTON_CONTROLLER_HH
#define DUMUX_NONEQUILIBRIUM_NEWTON_CONTROLLER_HH

#include <algorithm>

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief A nonequilibrium specific controller for the newton solver.
 * This controller calls the velocity averaging in the problem after each iteration.
 */
template <class Scalar,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class NonEquilibriumNewtonController : public NewtonController<Scalar, Comm>
{
    using ParentType = NewtonController<Scalar, Comm>;

public:
    using ParentType::ParentType;

    template<class JacobianAssembler, class SolutionVector>
    void newtonUpdate(JacobianAssembler& assembler,
                      SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        ParentType::newtonUpdate(assembler,
                                 uCurrentIter,
                                 uLastIter,
                                 deltaU);

        auto& gridVariables = assembler.gridVariables();
        // Averages the face velocities of a vertex. Implemented in the model.
        // The velocities are stored in the model.
        gridVariables.calcVelocityAverage(uCurrentIter);
    }
};

} // end namespace Dumux
#endif // DUMUX_VELO_PROB_AVERAGE_NEWTON_CONTROLLER_HH
