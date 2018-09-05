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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVariables
 */
#ifndef DUMUX_FREEFLOW_RANS_ZEROEQ_GRID_VARIABLES_HH
#define DUMUX_FREEFLOW_RANS_ZEROEQ_GRID_VARIABLES_HH

#include <dumux/freeflow/rans/gridvariables.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class storing data associated to scvs and scvfs
 * \tparam GG the type of the grid geometry
 * \tparam GVV the type of the grid volume variables
 * \tparam GFVC the type of the grid flux variables cache
 * \tparam GFV the type of the grid face variables
 */
template<class NavierStokesGridVariables>
class ZeroEqGridVariables : public RANSGridVariables<NavierStokesGridVariables, ZeroEqGridVariables<NavierStokesGridVariables>>
{
    using ThisType = ZeroEqGridVariables<NavierStokesGridVariables>;
    using ParentType = RANSGridVariables<NavierStokesGridVariables, ThisType>;
    using FVGridGeometry = typename ParentType::GridGeometry;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto dim = FVGridGeometry::GridView::dimension;

    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<typename ParentType::Scalar, dim, dim>;

    // static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    // static constexpr auto faceIdx = FVGridGeometry::faceIdx();

public:
    //! Export the Scalar type
    using Scalar = typename ParentType::Scalar;

    using ParentType::ParentType;





};

} // end namespace Dumux

#endif
