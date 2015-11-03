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
 * \brief A MpNc specific controller for the newton solver.
 *
 *        This controller calls the velocity averaging in the model after each iteration.
 */
#ifndef DUMUX_VELO_MODEL_NEWTON_CONTROLLER_HH
#define DUMUX_VELO_MODEL_NEWTON_CONTROLLER_HH

#include <algorithm>

#include <dumux/nonlinear/newtoncontroller.hh>
#include "mpncnewtoncontroller.hh"
#include "mpncproperties.hh"

namespace Dumux {
/*!
 * \brief A kinetic-MpNc specific controller for the newton solver.
 *
 * This controller calls the velocity averaging in the problem after each iteration.
 *
 * Everything else is taken from the standard mpnc newtoncontroller.
 */
template <class TypeTag>
class VeloModelNewtonController : public MPNCNewtonController<TypeTag>
{
    typedef MPNCNewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum {velocityAveragingInModel  = GET_PROP_VALUE(TypeTag, VelocityAveragingInModel)};

public:
    VeloModelNewtonController(const Problem &problem)
        : ParentType(problem)
    {};

    void newtonBeginStep(){
        ParentType::newtonBeginStep();

        // Averages the face velocities of a vertex. Implemented in the model.
        // The velocities are stored in the model.
        if(velocityAveragingInModel)
            this->problem_().model().calcVelocityAverage();
    }

    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        ParentType::newtonUpdate(uCurrentIter,
                                uLastIter,
                                deltaU);
    }
};

} // end namespace Dumux
#endif // DUMUX_VELO_PROB_AVERAGE_NEWTON_CONTROLLER_HH
