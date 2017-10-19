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
 * \brief A one-phase n-component specific controller for the newton solver.
 */
#ifndef DUMUX_1PNC_NEWTON_CONTROLLER_HH
#define DUMUX_1PNC_NEWTON_CONTROLLER_HH

#include "properties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup Newton
 * \ingroup OnePNCModel
 * \brief A one-phase n-component specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class OnePNCNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;


public:
  OnePNCNewtonController(Problem &problem)
      : ParentType(problem)
    {};

    /*!
     * \brief
     * Suggest a new time step size based either on the number of newton
     * iterations required or on the variable switch
     *
     * \param uCurrentIter The current global solution vector
     * \param uLastIter The previous global solution vector
     *
     */
//     void newtonEndStep(SolutionVector &uCurrentIter,
//                        const SolutionVector &uLastIter)
//     {
//     }

    bool newtonConverged()
    {

        return ParentType::newtonConverged();
    }
};
}

#endif
