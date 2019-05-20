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
 * \ingroup Discretization
 * \brief Base class for all finite volume grid geometries
 */
#ifndef DUMUX_DISCRETIZATION_BASE_FV_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BASE_FV_GRID_GEOMETRY_HH

#include "basegridgeometry.hh"

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Base class for all finite volume grid geometries
 * \tparam Impl the type of the actual implementation
 * \tparam GV the grid view type
 * \tparam Traits the fv geometry traits
 * \note This has been deprecated in favour of a discretization-scheme
 *       agnostic implementation in Dumux::BaseGridGeometry.
 */
template<class Impl, class GV, class Traits>
class DUNE_DEPRECATED_MSG("Use BaseGridGeometry instead. Will be removed after 3.1!")
BaseFVGridGeometry : public BaseGridGeometry<GV, Traits>
{
    using ParentType = BaseGridGeometry<GV, Traits>;
    using Element = typename GV::template Codim<0>::Entity;

public:
    BaseFVGridGeometry(const GV& gridView)
    : BaseGridGeometry<GV, Traits>(gridView) {}
};

} // end namespace Dumux

#endif
