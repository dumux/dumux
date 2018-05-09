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
 * \brief Base class for all geomechanical problems
 */
#ifndef DUMUX_GEOMECHANICS_FV_PROBLEM_HH
#define DUMUX_GEOMECHANICS_FV_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup Geomechanics
 * \brief Base class for all geomechanical problems
 * \note We actually require the same functionality as the
 *       porous medium flow problem, which is why we simply
 *       inherit from that here.
 */
template<class TypeTag>
class GeomechanicsFVProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

public:
    //! pull up the constructor of the parent class
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
