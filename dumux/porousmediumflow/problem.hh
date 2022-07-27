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
 * \ingroup PorousmediumflowModels
 * \brief Base class for all porous media problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_PROBLEM_HH
#define DUMUX_POROUS_MEDIUM_FLOW_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Base class for all fully implicit porous media problems.
 *
 * TODO: derive from base problem property?
 */
template<class TypeTag>
class PorousMediumFlowProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:

    //! Use constructors of the base class
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
