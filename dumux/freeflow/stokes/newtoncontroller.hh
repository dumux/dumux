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
 * \brief A Stokes specific controller for the newton solver.
 */
#ifndef DUMUX_STOKES_NEWTON_CONTROLLER_HH
#define DUMUX_STOKES_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup BoxStokesModel
 * \ingroup Newton
 * \brief A Stokes-specific controller for the nonlinear Newton solver, which sets
 *        different parameters for the relative tolerance, target steps and maximum steps.
 */
template <class TypeTag>
class StokesNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    StokesNewtonController(const Problem &problem)
        : ParentType(problem)
    {
        Dune::FMatrixPrecision<>::set_singular_limit(1e-35);
    }
};
} // end namespace Dumux

#endif
