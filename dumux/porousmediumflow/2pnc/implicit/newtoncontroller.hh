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
 * \brief A two-phase n-component specific controller for the newton solver.
 */
#ifndef DUMUX_2PNC_NEWTON_CONTROLLER_HH
#define DUMUX_2PNC_NEWTON_CONTROLLER_HH

#include "model.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup Newton
 * \ingroup TwoPNCModel
 * \brief A two-phase n-component specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class TwoPNCNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

public:
    TwoPNCNewtonController(Problem &problem) : ParentType(problem) {};

    /*!
     * \brief Returns true if the current solution can be considered to
     *        be accurate enough
     */
    bool newtonConverged()
    {
        if (this->method().model().switched())
            return false;

        return ParentType::newtonConverged();
    }
};

} // end namespace Dumux

#endif
